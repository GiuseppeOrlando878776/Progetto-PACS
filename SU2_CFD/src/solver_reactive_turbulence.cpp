#include "../include/solver_reactive_turbulence.hpp"

#include "../../Common/include/su2_assert.hpp"
#include "../include/variable_reactive.hpp"

//
//
/*!
 * \brief Turbulence model class constructor
 */
//
//
CReactiveTurbSolver::CReactiveTurbSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh): CTurbSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;

  bool restart = (config->GetRestart() || config->GetRestart_Flow());

  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  nVar = 2;
  nSpecies = CReactiveEulerVariable::GetnSpecies();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/
  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];

  /*--- Single grid simulation ---*/
  if(iMesh == MESH_0) {
    /*--- Define some auxiliary vector related with the residual ---*/
    Residual     = new su2double[nVar];
    Residual_RMS = new su2double[nVar];
    Residual_i   = new su2double[nVar];
    Residual_j   = new su2double[nVar];
    Residual_Max = new su2double[nVar];

    /*--- Define some structures for locating max residuals ---*/
    Point_Max = new unsigned long[nVar];
    Point_Max_Coord = new su2double*[nVar];
    for(iVar = 0; iVar < nVar; ++iVar)
      Point_Max_Coord[iVar] = new su2double[nDim];

    /*--- Define some auxiliary vector related with the solution ---*/
    Solution   = new su2double[nVar];
    Solution_i = new su2double[nVar];
    Solution_j = new su2double[nVar];

    /*--- Define some auxiliary vector related with the geometry ---*/
    Vector_i = new su2double[nDim];
    Vector_j = new su2double[nDim];

    /*--- Define some auxiliary vector related with the flow solution ---*/
    FlowPrimVar_i = new su2double [nSpecies + nDim + 5];
    FlowPrimVar_j = new su2double [nSpecies + nDim + 5];

    /*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for(iVar = 0; iVar < nVar; ++iVar) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/
    if(rank == MASTER_NODE)
      std::cout << "Initialize Jacobian structure (k-w model)." << std::endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if(config->GetKind_Linear_Solver_Prec() == LINELET ||
       config->GetKind_Linear_Solver() == SMOOTHER_LINELET) {
         nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
         if(rank == MASTER_NODE)
          std::cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << std::endl;
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }

  /*--- Computation of gradients by least squares ---*/
  if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for(iVar = 0; iVar < nVar; ++iVar)
    Cvector[iVar] = new su2double [nDim];
  }

  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];

  lowerlimit[CReactiveTurbVariable::TURB_KINE_INDEX] = 1.0e-10;
  upperlimit[CReactiveTurbVariable::TURB_KINE_INDEX] = 1.0e10;

  lowerlimit[CReactiveTurbVariable::TURB_OMEGA_INDEX] = 1.0e-4;
  upperlimit[CReactiveTurbVariable::TURB_OMEGA_INDEX] = 1.0e15;

  /*--- Flow infinity initialization stuff ---*/
  auto rhoInf    = config->GetDensity_FreeStreamND();
  auto muLamInf  = config->GetViscosity_FreeStreamND();
  auto Intensity = config->GetTurbulenceIntensity_FreeStream();
  auto viscRatio = config->GetTurb2LamViscRatio_FreeStream();
  auto VelInf    = config->GetVelocity_FreeStreamND();

  su2double VelMag = std::inner_product(VelInf, VelInf + nDim, VelInf, 0.0);
  VelMag = std::sqrt(VelMag);

  kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
  omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);

  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  su2double muT_Inf = rhoInf*kine_Inf/omega_Inf;

  /*--- Restart the solution from file information ---*/
  if(!restart || iMesh != MESH_0) {
    for(iPoint = 0; iPoint < nPoint; ++iPoint)
      node[iPoint] = new CReactiveTurbVariable(kine_Inf, omega_Inf, muT_Inf, nDim, nVar, config);
  }
  else
    LoadRestart(geometry, config, rank, muT_Inf);

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
}


//
//
/*!
 * \brief Load the file for restart
 */
//
//
void CReactiveTurbSolver::LoadRestart(CGeometry* geometry, CConfig* config, int rank, su2double muT_Inf) {
  unsigned long iPoint, index;
  su2double dull_val;
  unsigned short nZone = geometry->GetnZone();
  unsigned short iZone = config->GetiZone();
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);

  /*--- Restart the solution from file information ---*/
  std::ifstream restart_file;
  std::string filename = config->GetSolution_FlowFileName();

  /*--- Modify file name for multizone problems ---*/
  if(nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/
  if(dual_time || time_stepping) {
    int Unst_RestartIter;
    if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
    else if(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Open the restart file, throw an error if this fails. ---*/
  restart_file.open(filename.data(), std::ios::in);
  if(restart_file.fail()) {
    std::cout << "There is no turbulent restart file!!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*--- In case this is a parallel simulation, we need to perform the
        Global2Local index transformation first. ---*/
  std::map<unsigned long,unsigned long> Global2Local;

  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

  /*--- Read all lines in the restart file ---*/
  long iPoint_Local;
  unsigned long iPoint_Global;
  std::string text_line;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0,
  sbuf_NotMatching = 0;

  /*--- The first line is the header ---*/
  std::getline(restart_file, text_line);

  for(iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); ++iPoint_Global) {
    std::getline(restart_file, text_line);
    std::istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
          on the current processor, we will load and instantiate the vars. ---*/
    auto MI = Global2Local.find(iPoint_Global);
    if(MI != Global2Local.end()) {
      iPoint_Local = Global2Local[iPoint_Global];
      if(nDim == 2)
        point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >>
                      Solution[CReactiveTurbVariable::TURB_KINE_INDEX] >> Solution[CReactiveTurbVariable::TURB_OMEGA_INDEX];
      if(nDim == 3)
        point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >>
                      Solution[CReactiveTurbVariable::TURB_KINE_INDEX] >> Solution[CReactiveTurbVariable::TURB_OMEGA_INDEX];

      /*--- Instantiate the solution at this node, note that the muT_Inf should recomputed ---*/
      node[iPoint_Local] = new CReactiveTurbVariable(Solution[CReactiveTurbVariable::TURB_KINE_INDEX],
                                                     Solution[CReactiveTurbVariable::TURB_OMEGA_INDEX], muT_Inf, nDim, nVar, config);
      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/
  if(iPoint_Global_Local < nPointDomain)
    sbuf_NotMatching = 1;

  #ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
  #else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
  #endif
  if(rbuf_NotMatching != 0) {
    if(rank == MASTER_NODE) {
      std::cout << std::endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << std::endl;
      std::cout << "It could be empty lines at the end of the file." << std::endl << std::endl;
    }
    #ifndef HAVE_MPI
      std::exit(EXIT_FAILURE);
    #else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    #endif
  }

  /*--- Instantiate the variable class with an arbitrary solution
        at any halo/periodic nodes. The initial solution can be arbitrary,
        because a send/recv is performed immediately in the solver. ---*/
  for(iPoint = nPointDomain; iPoint < nPoint; ++iPoint)
    node[iPoint] = new CReactiveTurbVariable(Solution[CReactiveTurbVariable::TURB_KINE_INDEX],
                                             Solution[CReactiveTurbVariable::TURB_OMEGA_INDEX], muT_Inf, nDim, nVar, config);

  /*--- Close the restart file ---*/
  restart_file.close();
}


//
//
/*!
 * \brief Preprocessing
 */
//
//
void CReactiveTurbSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                        unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint;
  unsigned long ExtIter = config->GetExtIter();
  bool limiter_flow     = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter());

  for(iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);
  }

  /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
  if(config->GetKind_Gradient_Method() == GREEN_GAUSS)
    SetSolution_Gradient_GG(geometry, config);
  else if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    SetSolution_Gradient_LS(geometry, config);

  if(config->GetSpatialOrder() == SECOND_ORDER_LIMITER)
    SetSolution_Limiter(geometry, config);

  if(limiter_flow)
    solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
}

//
//
/*!
 * \brief Postprocessing
 */
//
//
void CReactiveTurbSolver::Postprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh) {
  su2double rho, omega, kine, muT;
  unsigned long iPoint;

  /*--- Compute mean flow and turbulence gradients ---*/
  if(config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    //solver_container[FLOW_SOL]->SetPrimitive_Gradient_GG(geometry, config);
    SetSolution_Gradient_GG(geometry, config);
  }
  else if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    //solver_container[FLOW_SOL]->SetPrimitive_Gradient_LS(geometry, config);
    SetSolution_Gradient_LS(geometry, config);
  }

  for(iPoint = 0; iPoint < nPoint; ++iPoint) {
    rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();

    /*--- Compute the eddy viscosity ---*/
    kine  = node[iPoint]->GetSolution(CReactiveTurbVariable::TURB_KINE_INDEX);
    omega = node[iPoint]->GetSolution(CReactiveTurbVariable::TURB_OMEGA_INDEX);
    muT = rho*kine/omega;
    node[iPoint]->SetmuT(muT);
  }
}

//
//
/*!
 * \brief Computing viscous residual
 */
//
//
void CReactiveTurbSolver::Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                           CConfig* config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;

  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Points coordinates, and normal vector ---*/
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Conservative variables w/o reconstruction ---*/
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive());
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    numerics->SetEddyViscosity(node[iPoint]->GetmuT(), node[jPoint]->GetmuT());

    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    numerics->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());

    /*--- Compute residual, and Jacobians ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Add and subtract residual, and update Jacobians ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
  }
}

//
//
/*!
 * \brief Computing source residual
 */
//
//
void CReactiveTurbSolver::Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                          CNumerics* second_numerics, CConfig* config, unsigned short iMesh) {
  for(unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Conservative variables w/o reconstruction ---*/
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);

    /*--- Gradient of the primitive and conservative variables ---*/
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

    /*--- Set volume ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Set distance to the surface ---*/
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

    /*--- Set vorticity and strain rate magnitude ---*/
    numerics->SetEddyViscosity(node[iPoint]->GetmuT(), 0.0);
    numerics->SetVorticity(solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(), NULL);
    numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(), 0.0);

    /*--- Compute the source term ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);

    /*--- Subtract residual and the Jacobian ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
  }
}

//
//
/*!
 * \brief HeatFlux Wall Boundary Condition
 */
//
//
void CReactiveTurbSolver::BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                           CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;

  su2double distance, density, laminar_viscosity;

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for(iDim = 0; iDim < nDim; ++iDim) {
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
                    (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = std::sqrt(distance);

      /*--- Set wall values ---*/
      density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
      laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();

      const su2double beta_1 = 3.0/40;
      Solution[CReactiveTurbVariable::TURB_KINE_INDEX] = 0.0;
      Solution[CReactiveTurbVariable::TURB_OMEGA_INDEX] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);

      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for(iVar = 0; iVar < nVar; ++iVar) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  } /*--- End of iVertex for loop ---*/
}

//
//
/*!
 * \brief Isothermal Wall Boundary Condition
 */
//
//
void CReactiveTurbSolver::BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                             CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

//
//
/*!
 * \brief Far Field Boundary Condition
 */
//
//
void CReactiveTurbSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                       CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;

  bool grid_movement = config->GetGrid_Movement();

  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Allocate the value at the infinity ---*/
      auto V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      auto V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/
      for(iVar = 0; iVar < nVar; ++iVar)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      Solution_j[CReactiveTurbVariable::TURB_KINE_INDEX] = kine_Inf;
      Solution_j[CReactiveTurbVariable::TURB_OMEGA_INDEX] = omega_Inf;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set Normal (it is necessary to change the sign) ---*/
      auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      std::transform(Normal, Normal + nDim, Normal, std::negate<su2double>());
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute residuals and Jacobians ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Add residuals and Jacobians ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  } /*--- End of iVertex for loop ---*/
}

//
//
/*!
 * \brief Inlet Boundary Condition
 */
//
//
void CReactiveTurbSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                   CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned short iVar;
  unsigned long iVertex, iPoint, Point_Normal;

  bool grid_movement  = config->GetGrid_Movement();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      std::transform(Normal, Normal + nDim, Normal, std::negate<su2double>());

      /*--- Allocate the value at the inlet ---*/
      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      auto V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states. Use free-stream l-w
            values for the turbulent state at the inflow. ---*/
      for(iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      Solution_j[CReactiveTurbVariable::TURB_KINE_INDEX]= kine_Inf;
      Solution_j[CReactiveTurbVariable::TURB_OMEGA_INDEX]= omega_Inf;

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    }
  }
}

//
//
/*!
 * \brief Outlet Boundary Condition
 */
//
//
void CReactiveTurbSolver::BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                    CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar;

  bool grid_movement  = config->GetGrid_Movement();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Allocate the value at the outlet ---*/
      auto V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      auto V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
            that the turbulent variable is copied from the interior of the
            domain to the outlet before computing the residual.
            Solution_i --> TurbVar_internal,
            Solution_j --> TurbVar_outlet ---*/
      for(iVar = 0; iVar < nVar; ++iVar) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set Normal (negate for outward convention) ---*/
      auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      std::transform(Normal, Normal + nDim, Normal, std::negate<su2double>());
      conv_numerics->SetNormal(Normal);

      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);

      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    }
  }
}

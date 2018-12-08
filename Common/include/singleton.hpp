#ifndef SU2_SINGLETON
#define SU2_SINGLETON

namespace Common {

  /*!
    * \brief This class provides a clear and clean way to specify, through inheritance, that a certain object can't be copied
    * \author G. Orlando
    */

  template<class T>
  class Singleton {
  public:
    /*!
      * \brief Default constructor
      */
    Singleton() = default;

    /*!
      * \brief Virtual destructor
      */
    virtual ~Singleton() = default;

    /*!
      * \brief Deleted copy constructor to prevent copies
      */
    Singleton(const Singleton&) = delete;

    /*!
      * \brief Deleted assignment operator to prevent copy-assignment
      */
    Singleton& operator=(const Singleton&) = delete;

  }; /*-- End of class Singleton ---*/

} /*-- End of Namespace Common ---*/

#endif

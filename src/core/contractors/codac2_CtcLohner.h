/** 
 *  \file codac2_CtcLohner.h
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Auguste Bourgois, Maël Godard
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_TimePropag.h"
#include "codac2_AnalyticFunction.h"
#include "codac2_TDomain.h"
#include "codac2_SlicedTube.h"

namespace codac2
{ 

  /**
   * \class GlobalEnclosureError
   *
   * \brief Encapsulates runtime error for global enclosure estimation failure
   */
  struct GlobalEnclosureError : public std::runtime_error {
    GlobalEnclosureError() : std::runtime_error(
        "Exceeded loop maximum range while looking for a global enclosure for the system.") {}
  };

  /**
   * \class LohnerAlgorithm
   *
   * \brief Simple first order Lohner algorithm to perform guaranteed integration of a system
   * \f$\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x})\f$
   */
  class LohnerAlgorithm {
  public:
    constexpr static const double FWD = 1., BWD = -1.;

    /**
     * \brief Creates a Lohner algorithm object
     *
     * \param f function defining the system \f$\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x})\f$
     * \param h time step of the integration method
     * \param forward forward integration if true, backaward otherwise
     * \param u0 initial condition of the system
     * \param contractions number of contractions of the global enclosure by the estimated local enclosure
     * \param eps inflation parameter for the global enclosure
     */
    LohnerAlgorithm(const AnalyticFunction<VectorType> *f,
                    double h,
                    bool forward,
                    const IntervalVector &u0 = IntervalVector::empty(1),
                    int contractions = 1,
                    double eps = 0.1);
    /**
     * \brief integrate the system over a given number of steps
     *
     * \param steps number of steps to integrate
     * \param H parameter to overwrite the integration time step
     * \return enclosure of the system's state
     */
    const IntervalVector &integrate(unsigned int steps, double H = -1);

    /**
     * \brief contract the global & local enclosure of the previous integration step
     *
     * \param x estimation of the global enclosure
     */
    void contractStep(const IntervalVector &x);

    /**
     * \brief Returns the current global enclosure, i.e. the box enclosing the trajectories between times \f$k-1\f$ and
     * \f$k\f$
     *
     * \return global enclosure
     */
    const IntervalVector &getLocalEnclosure() const;

    /**
     * \brief Returns the current local enclosure, i.e. the box enclosing the trajectories at time \f$k\f$
     *
     * \return local enclosure \f$[\mathbf{x}_{k}]\f$
     */
    const IntervalVector &getGlobalEnclosure() const;
  private:
    /**
     * \brief Computes an estimation of the global enclosure of the system
     *
     * \param initialGuess initial guess for the global enclosure
     * \param dir forward or backward in time
     * \return estimation of the global enclosure
     */
    IntervalVector globalEnclosure(const IntervalVector &initialGuess, double dir);

    unsigned int dim; //!< dimension of the system's state
    double h; //!< integration time step
    double direction; //!< forward or backward integration
    double eps; //!< inflation parameter for the global enclosure
    int contractions; //!< number of contractions of the global enclosure by the estimated local enclosure
    IntervalVector u; //!< local enclosure
    IntervalVector z; //!< Taylor-Lagrange remainder (order 2)
    IntervalVector r; //!< enclosure of uncertainties in the frame given by the matrix B
    IntervalVector u_tilde; //!< global enclosure
    Matrix B, Binv;
    Vector u_hat; //!< center of the box u
    const AnalyticFunction<VectorType> *f; //!< litteral function of the system \f$\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x})\f$
  };

  class CtcLohner
  {
    public:

      CtcLohner(const AnalyticFunction<VectorType>& f, int contractions = 5, double eps = 0.1);

      void contract(SlicedTube<IntervalVector> &tube, TimePropag t_propa = TimePropag::FWD_BWD) const;


    protected:

      AnalyticFunction<VectorType> m_f; //!< forward function
      int contractions; //!< number of contractions of the global enclosure by the estimated local enclosure
      int dim; //!< dimension of the state vector
      double eps; //!< inflation parameter for the global enclosure

      static const std::string m_ctc_name; //!< class name (mainly used for CN Exceptions)
      static std::vector<std::string> m_str_expected_doms;
  };
}
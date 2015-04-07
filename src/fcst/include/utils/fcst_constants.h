// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: fcst_constants.h
// - Description: This namespace contains all necessary FCST constants
// - Developers: Valentin N. Zingan,    University of Alberta
//               Marc Secanell Gallart, University of Alberta
// - Id: $Id: fcst_constants.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FCST_CONSTANTS_H_
#define _FCST_CONSTANTS_H_

#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>

const unsigned int dim = deal_II_dimension;

namespace Constants
{
      // --- GROUP 1 ---

      /**
       * Universal gas constant, \f$ R = 8.314462176 \quad \left[ \frac{\text{J}}{\text{mol K}} \right] \f$.
       */
      extern inline double R() { return 8.314462176; }

      /**
       * Faraday constant, \f$ F = 9.648533992 \cdot 10^4 \quad \left[ \frac{\text{C}}{\text{mol}} \right] \f$.
       */
      extern inline double F() { return 9.648533992e4; }

      /**
       * \f$ \pi = 3.141592654 \f$.
       */
      extern inline double Pi() { return 3.141592654; }

      /**
       * Permittivity of free space, \f$ \epsilon_0 = 8.854187818 \cdot 10^{-12} \quad \left[ \frac{\text{F}}{\text{m}} \right] \f$.
       */
      extern inline double E0() { return 8.854187818e-12; }

      /**
       * Boltzmann constant, \f$ k = 8.617332478 \cdot 10^{-5} \quad \left[ \frac{\text{eV}}{\text{K}} \right] \f$.
       */
      extern inline double K() { return 8.617332478e-5; }

      // --- GROUP 2 ---

      /**
       * Gravitational acceleration, \f$ \mathbf{g} = \{ g_{\alpha} \}_{\alpha = 1}^d \quad \text{such that} \quad \forall \alpha \neq d : \quad g_{\alpha} = 0 \quad \left[ \frac{\text{m}}{\text{sec}^2} \right] \quad
       * \text{and} \quad g_d = -9.81 \quad \left[ \frac{\text{m}}{\text{sec}^2} \right] \f$.
       */
      extern inline const dealii::Tensor<1,dim> gravity_acceleration()
      {
             dealii::Tensor<1,dim> result;
             result[dim-1] = -9.81;
             return result;
      }

      /**
       * Unit tensor, \f$ \hat{ \mathbf{I} } = \{ \delta_{\alpha \beta} \}_{\alpha,\beta = 1}^d \f$.
       */
      extern inline const dealii::SymmetricTensor<2,dim> unit_tensor()
      {
             dealii::SymmetricTensor<2,dim> result;

             for(unsigned int i = 0; i < dim; ++i)
                    for(unsigned int j = 0; j < dim; ++j)
                           if( i == j )
                                  result[i][j] = 1.0;

             return result;
      }

      // --- GROUP 3 ---

      /**
       * Coefficient \f$ A = 1.16145 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (viscosity and thermal conductivity).
       */
      extern inline double A_vk() { return 1.16145; }

      /**
       * Coefficient \f$ B = 0.14874 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (viscosity and thermal conductivity).
       */
      extern inline double B_vk() { return 0.14874; }

      /**
       * Coefficient \f$ C = 0.52487 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (viscosity and thermal conductivity).
       */
      extern inline double C_vk() { return 0.52487; }

      /**
       * Coefficient \f$ D = 0.77320 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (viscosity and thermal conductivity).
       */
      extern inline double D_vk() { return 0.77320; }

      /**
       * Coefficient \f$ E = 2.16178 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (viscosity and thermal conductivity).
       */
      extern inline double E_vk() { return 2.16178; }

      /**
       * Coefficient \f$ F = 2.43787 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (viscosity and thermal conductivity).
       */
      extern inline double F_vk() { return 2.43787; }

      // --- GROUP 4 ---

      /**
       * Coefficient \f$ A = 1.06036 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double A_diff() { return 1.06036; }

      /**
       * Coefficient \f$ B = 0.15610 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double B_diff() { return 0.15610; }

      /**
       * Coefficient \f$ C = 0.19300 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double C_diff() { return 0.19300; }

      /**
       * Coefficient \f$ D = 0.47635 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double D_diff() { return 0.47635; }

      /**
       * Coefficient \f$ E = 1.03587 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double E_diff() { return 1.03587; }

      /**
       * Coefficient \f$ F = 1.52996 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double F_diff() { return 1.52996; }

      /**
       * Coefficient \f$ G = 1.76474 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double G_diff() { return 1.76474; }

      /**
       * Coefficient \f$ H = 3.89411 \f$ of the Neufeld-Jansen-Aziz
       * collision integrals formula (diffusivity).
       */
      extern inline double H_diff() { return 3.89411; }

      // --- GROUP 5 ---

      /**
       * Coefficient \f$ b_0 = -2.1794               \f$ of the Springer-Zawodzinski-Gottesfeld
       * water vapor saturation pressure formula.
       */
      extern inline double b_0() { return -2.1794;    }

      /**
       * Coefficient \f$ b_1 = 2.9530 \cdot 10^{-2}  \f$ of the Springer-Zawodzinski-Gottesfeld
       * water vapor saturation pressure formula.
       */
      extern inline double b_1() { return 2.9530e-2;  }

      /**
       * Coefficient \f$ b_2 = -9.1837 \cdot 10^{-5} \f$ of the Springer-Zawodzinski-Gottesfeld
       * water vapor saturation pressure formula.
       */
      extern inline double b_2() { return -9.1837e-5; }

      /**
       * Coefficient \f$ b_3 = 1.4454 \cdot 10^{-7}  \f$ of the Springer-Zawodzinski-Gottesfeld
       * water vapor saturation pressure formula.
       */
      extern inline double b_3() { return 1.4454e-7;  }
}

#endif
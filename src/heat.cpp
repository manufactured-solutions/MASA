// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// heat.cpp: These are the MASA class member functions and constructors
//          For the Heat Equation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *   steady -- constant coeff
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::heateq_1d_steady_const<Scalar>::heateq_1d_steady_const()
{
    this->mmsname = "heateq_1d_steady_const";
    this->dimension=1;

    this->register_var("A_x",&A_x);   
    this->register_var("k_0",&k_0);

}//done with constructor

template <typename Scalar>
int MASA::heateq_1d_steady_const<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("A_x",1.4);
  err += this->set_var("k_0",.82);
  
  return err;

}

template <typename Scalar>
Scalar MASA::heateq_1d_steady_const<Scalar>::eval_q_t(Scalar x)
{
  Scalar Q_T;
  Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

template <typename Scalar>
Scalar MASA::heateq_1d_steady_const<Scalar>::eval_exact_t(Scalar x)
{
  Scalar exact_t;
  exact_t = cos(A_x * x);
  return exact_t;
}

template <typename Scalar>
MASA::heateq_2d_steady_const<Scalar>::heateq_2d_steady_const()
{
    this->mmsname = "heateq_2d_steady_const";
    this->dimension=2;

    this->register_var("A_x",&A_x);   
    this->register_var("k_0",&k_0);
    this->register_var("B_y",&B_y);
    
}//done with constructor

template <typename Scalar>
int MASA::heateq_2d_steady_const<Scalar>::init_var()
{
  int err = 0;

  // default
  Scalar param=1.2;

  err += this->set_var("A_x",param);
  err += this->set_var("B_y",param);
  err += this->set_var("k_0",param);
  
  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_2d_steady_const<Scalar>::eval_q_t(Scalar x,Scalar y)
{
  Scalar Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

template <typename Scalar>
Scalar MASA::heateq_2d_steady_const<Scalar>::eval_exact_t(Scalar x,Scalar y)
{
  Scalar exact_t;
  exact_t = cos(A_x * x) * cos(B_y * y);
  return exact_t;
}

template <typename Scalar>
MASA::heateq_3d_steady_const<Scalar>::heateq_3d_steady_const()
{
    this->mmsname = "heateq_3d_steady_const";
    this->dimension=3;

    this->register_var("A_x",&A_x);   
    this->register_var("k_0",&k_0);
    this->register_var("B_y",&B_y);
    this->register_var("C_z",&C_z);

}//done with constructor

template <typename Scalar>
int MASA::heateq_3d_steady_const<Scalar>::init_var()
{
  int err = 0;
  
  // default
  Scalar param=1.2;


  err += this->set_var("A_x",param);
  err += this->set_var("B_y",param);
  err += this->set_var("k_0",param);
  err += this->set_var("C_z",param);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_3d_steady_const<Scalar>::eval_q_t(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  return Q_T;
}

template <typename Scalar>
Scalar MASA::heateq_3d_steady_const<Scalar>::eval_exact_t(Scalar x,Scalar y,Scalar z)
{
  Scalar exact_t;
  exact_t = cos(A_x * x) * cos(B_y * y) * cos(C_z * z);
  return exact_t;
}



/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *   Unsteady -- Const Coeff
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::heateq_1d_unsteady_const<Scalar>::heateq_1d_unsteady_const()
{
  this->mmsname = "heateq_1d_unsteady_const";
  this->dimension=1;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("D_t",&D_t);
  this->register_var("cp_0",&cp_0);
  this->register_var("A_t",&A_t);
  this->register_var("rho",&rho);

}//done with constructor

template <typename Scalar>
int MASA::heateq_1d_unsteady_const<Scalar>::init_var()
{
  int err=0;

  err += this->set_var("A_x",1.817);
  err += this->set_var("k_0",.1984);
  err += this->set_var("D_t",181.4);
  err += this->set_var("cp_0",.104);
  err += this->set_var("A_t", 12.4);
  err += this->set_var("rho",.2380);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_1d_unsteady_const<Scalar>::eval_q_t(Scalar x,Scalar t)
{
  Scalar Q_T = cos(A_x * x + A_t * t) * cos(D_t * t) * k_0 * A_x * A_x - (sin(A_x * x + A_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(D_t * t) * D_t) * rho * cp_0;
  return Q_T;
}


template <typename Scalar>
MASA::heateq_2d_unsteady_const<Scalar>::heateq_2d_unsteady_const()
{
  this->mmsname = "heateq_2d_unsteady_const";
  this->dimension=2;
    
  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("D_t",&D_t);
  this->register_var("cp_0",&cp_0);
  this->register_var("A_t",&A_t);
  this->register_var("rho",&rho);
  this->register_var("B_y",&B_y);
  this->register_var("B_t",&B_t);

}//done with constructor

template <typename Scalar>
int MASA::heateq_2d_unsteady_const<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("A_x",1.1);   
  err += this->set_var("k_0",5.3);
  err += this->set_var("D_t",7.8);
  err += this->set_var("cp_0",0.1);
  err += this->set_var("A_t",12.3);
  err += this->set_var("rho",1.2);
  err += this->set_var("B_y",5.4);
  err += this->set_var("B_t",8.09);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_2d_unsteady_const<Scalar>::eval_q_t(Scalar x,Scalar y, Scalar t)
{
  Scalar Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}

template <typename Scalar>
MASA::heateq_3d_unsteady_const<Scalar>::heateq_3d_unsteady_const()
{
  this->mmsname = "heateq_3d_unsteady_const";
  this->dimension=3;
  
  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("D_t",&D_t);
  this->register_var("cp_0",&cp_0);
  this->register_var("A_t",&A_t);
  this->register_var("rho",&rho);
  this->register_var("B_y",&B_y);
  this->register_var("B_t",&B_t);
  this->register_var("C_z",&C_z);
  this->register_var("C_t",&C_t);

}//done with constructor

template <typename Scalar>
int MASA::heateq_3d_unsteady_const<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("A_x",1.4);
  err += this->set_var("k_0",.82);
  err += this->set_var("D_t",.12);
  err += this->set_var("cp_0",3.1);
  err += this->set_var("A_t",0.01);
  err += this->set_var("rho",4.0);
  err += this->set_var("B_y",11.01);
  err += this->set_var("B_t",1.01);
  err += this->set_var("C_z",0.90);
  err += this->set_var("C_t",12.34);
  
  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_3d_unsteady_const<Scalar>::eval_q_t(Scalar x,Scalar y, Scalar z,Scalar t)
{
  Scalar Q_T = -(sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * A_t + cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * B_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * C_t + cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * D_t) * rho * cp_0 + (A_x * A_x + B_y * B_y + C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * k_0;
  return Q_T;
}

/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *   Unsteady -- Variable Coef
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::heateq_1d_unsteady_var<Scalar>::heateq_1d_unsteady_var()
{
  this->mmsname = "heateq_1d_unsteady_var";
  this->dimension=1;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("D_t",&D_t);
  this->register_var("cp_0",&cp_0);
  this->register_var("A_t",&A_t);
  this->register_var("rho",&rho);
  this->register_var("cp_1",&cp_1);
  this->register_var("cp_2",&cp_2);
  this->register_var("k_1",&k_1);
  this->register_var("k_2",&k_2);

}//done with constructor

template <typename Scalar>
int MASA::heateq_1d_unsteady_var<Scalar>::init_var()
{
  int err = 0;

  Scalar param =1.093;

  err += this->set_var("A_x",param);   
  err += this->set_var("k_0",param);
  err += this->set_var("D_t",param);
  err += this->set_var("cp_0",param);
  err += this->set_var("A_t",param);
  err += this->set_var("rho",param);
  err += this->set_var("cp_1",param);
  err += this->set_var("cp_2",param);
  err += this->set_var("k_1",param);
  err += this->set_var("k_2",param);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_1d_unsteady_var<Scalar>::eval_q_t(Scalar x, Scalar t)
{
  Scalar Q_T = -pow(cos(D_t * t), Scalar(2)) * A_x * A_x * k_1 + Scalar(3) * pow(cos(A_x * x + A_t * t), Scalar(3)) * pow(cos(D_t * t), Scalar(3)) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_0 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_0 + (A_x * A_x * k_0 - Scalar(2) * pow(cos(D_t * t), Scalar(2)) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_1 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_1) * cos(A_x * x + A_t * t) * cos(D_t * t) + (Scalar(2) * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(D_t * t) * A_t * rho * cp_2 - cos(A_x * x + A_t * t) * sin(D_t * t) * D_t * rho * cp_2) * pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2));
  return Q_T;
}

template <typename Scalar>
MASA::heateq_2d_unsteady_var<Scalar>::heateq_2d_unsteady_var()
{
  this->mmsname = "heateq_2d_unsteady_var";
  this->dimension=2;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("D_t",&D_t);
  this->register_var("cp_0",&cp_0);
  this->register_var("A_t",&A_t);
  this->register_var("rho",&rho);
  this->register_var("cp_1",&cp_1);
  this->register_var("cp_2",&cp_2);
  this->register_var("k_1",&k_1);
  this->register_var("k_2",&k_2);
  this->register_var("B_y",&B_y);
  this->register_var("B_t",&B_t);

}//done with constructor

template <typename Scalar>
int MASA::heateq_2d_unsteady_var<Scalar>::init_var()
{
  int err = 0;
  Scalar param =1.093;

  err += this->set_var("A_x",param);   
  err += this->set_var("k_0",param);
  err += this->set_var("D_t",param);
  err += this->set_var("cp_0",param);
  err += this->set_var("A_t",param);
  err += this->set_var("rho",param);
  err += this->set_var("cp_1",param);
  err += this->set_var("cp_2",param);
  err += this->set_var("k_1",param);
  err += this->set_var("k_2",param);
  err += this->set_var("B_y",param);
  err += this->set_var("B_t",param);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_2d_unsteady_var<Scalar>::eval_q_t(Scalar x,Scalar y,Scalar t)
{
  Scalar Q_T = -pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t * rho * cp_0 - pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * B_y * B_y * k_1 - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t * rho * cp_0 - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t * rho * cp_0 + Scalar(3) * pow(cos(A_x * x + A_t * t), Scalar(3)) * pow(cos(B_y * y + B_t * t), Scalar(3)) * pow(cos(D_t * t), Scalar(3)) * (A_x * A_x + B_y * B_y) * k_2 + (Scalar(2) * A_x * A_x * k_1 - sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t * rho * cp_2 + Scalar(2) * B_y * B_y * k_1 - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t * rho * cp_2 - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t * rho * cp_2) * pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) + (A_x * A_x * k_0 - Scalar(2) * pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * A_x * A_x * k_2 - sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t) * A_t * rho * cp_1 + B_y * B_y * k_0 - Scalar(2) * pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * B_y * B_y * k_2 - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(D_t * t) * B_t * rho * cp_1 - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(D_t * t) * D_t * rho * cp_1) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(D_t * t);
  return Q_T;
}

template <typename Scalar>
MASA::heateq_3d_unsteady_var<Scalar>::heateq_3d_unsteady_var()
{
  this->mmsname = "heateq_3d_unsteady_var";
  this->dimension=3;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("D_t",&D_t);
  this->register_var("cp_0",&cp_0);
  this->register_var("A_t",&A_t);
  this->register_var("rho",&rho);
  this->register_var("cp_1",&cp_1);
  this->register_var("cp_2",&cp_2);
  this->register_var("k_1",&k_1);
  this->register_var("k_2",&k_2);
  this->register_var("B_y",&B_y);
  this->register_var("B_t",&B_t);
  this->register_var("C_z",&C_z);
  this->register_var("C_t",&C_t);

}//done with constructor

template <typename Scalar>
int MASA::heateq_3d_unsteady_var<Scalar>::init_var()
{
  int err = 0;
  Scalar param =1.093;

  err += this->set_var("A_x",param);   
  err += this->set_var("k_0",param);
  err += this->set_var("D_t",param);
  err += this->set_var("cp_0",param);
  err += this->set_var("A_t",param);
  err += this->set_var("rho",param);
  err += this->set_var("cp_1",param);
  err += this->set_var("cp_2",param);
  err += this->set_var("k_1",param);
  err += this->set_var("k_2",param);
  err += this->set_var("B_y",param);
  err += this->set_var("B_t",param);
  err += this->set_var("C_z",param);
  err += this->set_var("C_t",param);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_3d_unsteady_var<Scalar>::eval_q_t(Scalar x,Scalar y,Scalar z,Scalar t)
{
  Scalar Q_T = -sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_0 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_0 * D_t - pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(C_z * z + C_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * k_1 * A_x * A_x - pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(C_z * z + C_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * k_1 * B_y * B_y - pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * k_1 * C_z * C_z + Scalar(3) * pow(cos(A_x * x + A_t * t), Scalar(3)) * pow(cos(B_y * y + B_t * t), Scalar(3)) * pow(cos(C_z * z + C_t * t), Scalar(3)) * pow(cos(D_t * t), Scalar(3)) * (A_x * A_x + B_y * B_y + C_z * C_z) * k_2 + (-sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_1 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_1 * D_t + k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - Scalar(2) * pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(C_z * z + C_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * k_2 * A_x * A_x - Scalar(2) * pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(C_z * z + C_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * k_2 * B_y * B_y - Scalar(2) * pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2)) * k_2 * C_z * C_z) * cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) + (-sin(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * A_t - cos(A_x * x + A_t * t) * sin(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * B_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * sin(C_z * z + C_t * t) * cos(D_t * t) * rho * cp_2 * C_t - cos(A_x * x + A_t * t) * cos(B_y * y + B_t * t) * cos(C_z * z + C_t * t) * sin(D_t * t) * rho * cp_2 * D_t + Scalar(2) * k_1 * A_x * A_x + Scalar(2) * k_1 * B_y * B_y + Scalar(2) * k_1 * C_z * C_z) * pow(cos(A_x * x + A_t * t), Scalar(2)) * pow(cos(B_y * y + B_t * t), Scalar(2)) * pow(cos(C_z * z + C_t * t), Scalar(2)) * pow(cos(D_t * t), Scalar(2));
  return Q_T;
}



/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *        steady
 *
 *        variable coefficients
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::heateq_1d_steady_var<Scalar>::heateq_1d_steady_var()
{
  this->mmsname = "heateq_1d_steady_var";
  this->dimension=1;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("k_1",&k_1);
  this->register_var("k_2",&k_2);

}//done with constructor

template <typename Scalar>
int MASA::heateq_1d_steady_var<Scalar>::init_var()
{
  int err = 0;
  Scalar param = 1.4;
  err += this->set_var("A_x",param);   
  err += this->set_var("k_0",param);
  err += this->set_var("k_1",param);
  err += this->set_var("k_2",param);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_1d_steady_var<Scalar>::eval_q_t(Scalar x)
{
  Scalar Q_T = Scalar(3) * A_x * A_x * k_2 * pow(cos(A_x * x), Scalar(3)) + Scalar(2) * A_x * A_x * k_1 * pow(cos(A_x * x), Scalar(2)) - A_x * A_x * k_1 + (k_0 - Scalar(2) * k_2) * A_x * A_x * cos(A_x * x);
  return Q_T;
}

template <typename Scalar>
MASA::heateq_2d_steady_var<Scalar>::heateq_2d_steady_var()
{
  this->mmsname = "heateq_2d_steady_var";
  this->dimension=2;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("k_1",&k_1);
  this->register_var("k_2",&k_2);
  this->register_var("B_y",&B_y);

}//done with constructor

template <typename Scalar>
int MASA::heateq_2d_steady_var<Scalar>::init_var()
{
  int err = 0;
  Scalar param = 1.4;
  err += this->set_var("A_x",param);   
  err += this->set_var("k_0",param);
  err += this->set_var("k_1",param);
  err += this->set_var("k_2",param);
  err += this->set_var("B_y",param);
  return err;

} // done with variable initializer


template <typename Scalar>
Scalar MASA::heateq_2d_steady_var<Scalar>::eval_q_t(Scalar x,Scalar y)
{
  Scalar Q_T = (Scalar(3) * A_x * A_x + Scalar(3) * B_y * B_y) * k_2 * pow(cos(A_x * x), Scalar(3)) * pow(cos(B_y * y), Scalar(3)) + (Scalar(2) * A_x * A_x + Scalar(2) * B_y * B_y) * k_1 * pow(cos(A_x * x), Scalar(2)) * pow(cos(B_y * y), Scalar(2)) - (pow(cos(B_y * y), Scalar(2)) * A_x * A_x + pow(cos(A_x * x), Scalar(2)) * B_y * B_y) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y - Scalar(2) * pow(cos(B_y * y), Scalar(2)) * k_2 * A_x * A_x - Scalar(2) * pow(cos(A_x * x), Scalar(2)) * k_2 * B_y * B_y) * cos(A_x * x) * cos(B_y * y);
  return Q_T;
}


template <typename Scalar>
MASA::heateq_3d_steady_var<Scalar>::heateq_3d_steady_var()
{
  this->mmsname = "heateq_3d_steady_var";
  this->dimension=3;

  this->register_var("A_x",&A_x);   
  this->register_var("k_0",&k_0);
  this->register_var("k_1",&k_1);
  this->register_var("k_2",&k_2);
  this->register_var("B_y",&B_y);
  this->register_var("C_z",&C_z);

}//done with constructor

template <typename Scalar>
int MASA::heateq_3d_steady_var<Scalar>::init_var()
{  
  int err = 0;
  Scalar param = 1.4;
  err += this->set_var("A_x",param);   
  err += this->set_var("k_0",param);
  err += this->set_var("k_1",param);
  err += this->set_var("k_2",param);
  err += this->set_var("B_y",param);
  err += this->set_var("C_z",param);
  return err;
} // done with variable initializer

template <typename Scalar>
Scalar MASA::heateq_3d_steady_var<Scalar>::eval_q_t(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_T = (Scalar(3) * A_x * A_x + Scalar(3) * B_y * B_y + Scalar(3) * C_z * C_z) * k_2 * pow(cos(A_x * x), Scalar(3)) * pow(cos(B_y * y), Scalar(3)) * pow(cos(C_z * z), Scalar(3)) + (Scalar(2) * A_x * A_x + Scalar(2) * B_y * B_y + Scalar(2) * C_z * C_z) * k_1 * pow(cos(A_x * x), Scalar(2)) * pow(cos(B_y * y), Scalar(2)) * pow(cos(C_z * z), Scalar(2)) - (pow(cos(B_y * y), Scalar(2)) * pow(cos(C_z * z), Scalar(2)) * A_x * A_x + pow(cos(A_x * x), Scalar(2)) * pow(cos(C_z * z), Scalar(2)) * B_y * B_y + pow(cos(A_x * x), Scalar(2)) * pow(cos(B_y * y), Scalar(2)) * C_z * C_z) * k_1 + (k_0 * A_x * A_x + k_0 * B_y * B_y + k_0 * C_z * C_z - Scalar(2) * pow(cos(B_y * y), Scalar(2)) * pow(cos(C_z * z), Scalar(2)) * k_2 * A_x * A_x - Scalar(2) * pow(cos(A_x * x), Scalar(2)) * pow(cos(C_z * z), Scalar(2)) * k_2 * B_y * B_y - Scalar(2) * pow(cos(A_x * x), Scalar(2)) * pow(cos(B_y * y), Scalar(2)) * k_2 * C_z * C_z) * cos(A_x * x) * cos(B_y * y) * cos(C_z * z);
  return Q_T;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::heateq_1d_steady_const);
MASA_INSTANTIATE_ALL(MASA::heateq_2d_steady_const);
MASA_INSTANTIATE_ALL(MASA::heateq_3d_steady_const);
MASA_INSTANTIATE_ALL(MASA::heateq_1d_steady_var);
MASA_INSTANTIATE_ALL(MASA::heateq_2d_steady_var);
MASA_INSTANTIATE_ALL(MASA::heateq_3d_steady_var);
MASA_INSTANTIATE_ALL(MASA::heateq_1d_unsteady_const);
MASA_INSTANTIATE_ALL(MASA::heateq_2d_unsteady_const);
MASA_INSTANTIATE_ALL(MASA::heateq_3d_unsteady_const);
MASA_INSTANTIATE_ALL(MASA::heateq_1d_unsteady_var);
MASA_INSTANTIATE_ALL(MASA::heateq_2d_unsteady_var);
MASA_INSTANTIATE_ALL(MASA::heateq_3d_unsteady_var);

// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
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
// $Author: nick $
// $Id: heat.cpp 18162 2011-03-01 05:23:07Z nick $
//
// radiation.cpp: These are the MASA class member functions and constructors
//          For Radiation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *   radiation_integrated_intensity
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::radiation_integrated_intensity<Scalar>::radiation_integrated_intensity()
{
    this->mmsname = "radiation_integrated_intensity";
    this->dimension=1;

    // registering variables
    this->register_var("no_gauss",&no_gauss);   

    // registering a vector
    this->register_vec("vec_mean",vec_mean);   
    this->register_vec("vec_amp",vec_amp);   
    this->register_vec("vec_stdev",vec_stdev);   
    this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::radiation_integrated_intensity<Scalar>::init_var()
{
  int err = 0;

  /*
  ! Random number generated :: Matlab
   a    =(/    7.4462631920000945d+01,    6.6885970491068861d+01,    &
    9.5002692320883099d+01,          6.8462339056010748d+01,         &
    5.5560137764689372d+01,          8.9012603416056891d+01,         &
    6.9486941848062671d+01,          6.2084564295691635d+01,         &
    7.0195607279405735d+01,          5.4822726258419429d+01,         &
    5.6598664630316755d+01,          9.7102529538774263d+01,         &
    9.7806727011490111d+01,          7.8760429753923276d+01,         &
    5.2988977147357794d+01,          6.1738995668620319d+01,         &
    6.7657928561103546d+01,          9.1059702009897961d+01,         &
    5.0770171882577756d+01,          5.2151190082890395d+01,         &
    5.8449501473135214d+01,          8.2455773747822604d+01,         &
    8.6586119282933510d+01,          8.2387298156815334d+01,         &
    7.2546185321547242d+01 /);

      b(1) = 0.2d0
      do i = 2, ngaussians
         b(i) = b(i-1)+(0.8d0-0.2d0)/(dble(ngaussians-1));
      enddo
    
      do i = 1, ngaussians
         c(i) = 0.05;
      enddo
      
  */


  // using marcos defaults

  err += this->set_var("no_gauss",25);

  // set size of vectors and set default values
  vec_mean.resize(no_gauss);
  vec_mean[0]=0.2;
  for(int it = 0;it<int(vec_amp.size());it++)
    {
      //vec_mean[it]=it;
      vec_mean[it]=vec_mean[it-1]+(0.8-0.2)/(24.0);;
    }
  
  vec_amp.resize(no_gauss);
  vec_amp[0] = 7.4462631920000945;
  vec_amp[1] = 6.6885970491068861;
  vec_amp[2] = 9.5002692320883099;
  vec_amp[3] = 6.8462339056010748;
  vec_amp[4] = 5.5560137764689372;
  vec_amp[5] = 8.9012603416056891;
  vec_amp[6] = 6.9486941848062671;
  vec_amp[7] = 6.2084564295691635;
  vec_amp[8] = 7.0195607279405735;
  vec_amp[9] = 5.4822726258419429;
  vec_amp[10]= 5.6598664630316755;
  vec_amp[11]= 9.7102529538774263;
  vec_amp[12]= 9.7806727011490111;
  vec_amp[13]= 7.8760429753923276;
  vec_amp[14]= 5.2988977147357794;
  vec_amp[15]= 6.1738995668620319;
  vec_amp[16]= 6.7657928561103546;
  vec_amp[17]= 9.1059702009897961;
  vec_amp[18]= 5.0770171882577756;
  vec_amp[19]= 5.2151190082890395;
  vec_amp[20]= 5.8449501473135214;
  vec_amp[21]= 8.2455773747822604;
  vec_amp[22]= 8.6586119282933510;
  vec_amp[23]= 8.2387298156815334;
  vec_amp[24]= 7.2546185321547242;

  /*
    for(int it = 0;it<int(vec_amp.size());it++)
    {
      //vec_amp[it]=12; //should sum to 60

    }
  */

  vec_stdev.resize(no_gauss);    
  vec_amp.resize(no_gauss);
  for(int it = 0;it<int(vec_amp.size());it++)
    {
      vec_stdev[it]=0.05; 
    }
 
  err += this->set_vec("vec_mean",vec_mean);
  err += this->set_vec("vec_amp",vec_amp);
  err += this->set_vec("vec_stdev",vec_stdev);

  return err;

}

template <typename Scalar>
int MASA::radiation_integrated_intensity<Scalar>::check_vec()
{

  if(vec_mean.size() != vec_amp.size() || vec_mean.size() != vec_stdev.size())
    {
      std::cout << "MASA WARNING:: in mms radiation_integrated_intensity--\n";
      std::cout << "Vectors vec_amp,vec_mean,vec_stdev are not identically sized!\n";	
      return 1;
    }

  return 0;
}

template <typename Scalar>
Scalar MASA::radiation_integrated_intensity<Scalar>::eval_q_u(Scalar x)
{
  // this is the manufactured solution: i.e. the gaussians contributions
  Scalar Q_I = 0;

  // error handling for vectors
  if(check_vec() == 1)
    {
      return -1;
    }
  
  // sum up intensity at particular location by
  // looping over gaussians of intensity
  for(int it = 0;it<int(vec_amp.size());it++)
    {
      // this is evaluating the gaussians contributions at 
      // a particular spatial location
      Q_I += vec_amp[it]*exp( -pow(x-vec_mean[it],2)/(2*pow(vec_stdev[it],2)));
    }

  return Q_I;  
}

template <typename Scalar>
Scalar MASA::radiation_integrated_intensity<Scalar>::eval_exact_u(Scalar x)
{
  // this is the source term: i.e. integrated intensity
  Scalar exact_I = 0;
  Scalar xs,s;
  // error handling for vectors
  if(check_vec() == 1)
    return -1;

  // error handling for x
  if(x <= 0)
    {
      std::cout << "MASA WARNING:: in mms radiation_integrated_intensity--\n";
      std::cout << "Integration bound must be greater than zero!\n";
      return -1;
    }

  // sum up gaussians for integrated intensity
  // achtung: some sort of iterator failure here, hacking together a loop for now
  //  for(std::vector<Scalar>::iterator it = vec_amp.begin(); it != vec_amp.end(); it++)
  for(int it = 0;it<int(vec_amp.size());it++)
    {
      xs = (x-vec_mean[it])/vec_stdev[it];
       s = -vec_mean[it]/vec_stdev[it];
      
       // integrate [0,inf]
       if(x>1000.0)
	 {
	   exact_I += vec_amp[it]*(1-phi(s));
	 }       
       else // integrate [0,x]
	 {
	   exact_I += vec_amp[it]*(phi(xs)-phi(s));
	 }
       
       // c = pow(2.0,0.5);
       // p = pow(pi,0.5);
       // 0.50*erf((c*vec_mean[it])/(2.0*vec_stdev[it]))*vec_amp[it]*sqrt(pi)*sqrt(2.0)*vec_stdev[it]  
       // - 0.5*erf((vec_stdev[it]*(-1.0+vec_mean[it]))/(2.0*vec_stdev[it])) * vec_amp[it]*p*c*vec_stdev[it];    
    }  
  
  return exact_I;
}

template <typename Scalar>
Scalar MASA::radiation_integrated_intensity<Scalar>::phi(Scalar x)
{
  Scalar out;
  // using c99 provided erf(double x) from cmath      
  out = 0.5 * (1 + erf(x));
  return out;  
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::radiation_integrated_intensity);

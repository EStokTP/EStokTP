C      Copyright (C) 2019 Carlo de Falco
C      
C      This program is free software: you can redistribute it and/or modify it
C      under the terms of the GNU General Public License as published by
C      the Free Software Foundation, either version 3 of the License, or
C      (at your option) any later version.
C      
C      This program is distributed in the hope that it will be useful, but
C      WITHOUT ANY WARRANTY; without even the implied warranty of
C      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C      GNU General Public License for more details.
C      
C      You should have received a copy of the GNU General Public License
C      along with this program.  If not, see <https://www.gnu.org/licenses/>.

C     THIS FILE DEFINES THE PROBLEM      
C     V1 = x^2 + y^2
C     V2 = x^2 + y^2 - (x + y -1)
C     MINIMIZE
C     V1 (x, y) = x^2 + y^2
C     CONSTRAINED BY
C     V1(x,y) - V2(x,y) = (x + y -1) = 0
C     EXPECTED CONSTRAINED MINIMU IS
C     x=y=1/2
      
      block data init
      implicit none
      integer n
      parameter (n = 2)
      double precision  qinit(n)
      common /general/ qinit
      data qinit / 0.0D0, 0.0D0/
      end
      
C     *****************************
      subroutine problem_size(res)
      implicit none
      integer n, res
      parameter (n = 2)
      
      res = n

      end subroutine

C     *****************************
      subroutine initial_guess(q)
      implicit none
      integer n
      parameter (n = 2)
      double precision q(n), qinit(n)
      common /general/ qinit
     
      q = qinit
      
      end subroutine
   

C     *****************************
      subroutine V1(q, z)
      implicit none
      integer n
      parameter (n = 2)
      double precision q(n), z
      
      z = q(1)**2 + q(2)**2

      end subroutine

C     *****************************
      subroutine V2(q, z)
      implicit none
      integer n
      parameter (n = 2)
      double precision q(n), z
      
      z = q(1)**2 + q(2)**2 - (q(1) + q(2) -1.0D0)

      end subroutine


C     *****************************
      subroutine dV1(q, dz)
      implicit none
      integer n
      parameter (n = 2)
      double precision q(n), dz(n)
      
      dz(1) = 2.0D0*q(1)
      dz(2) = 2.0D0*q(2)

      end subroutine

C     *****************************
      subroutine dV2(q, dz)
      implicit none
      integer n
      parameter (n = 2)
      double precision q(n), dz(n)
      
      dz(1) = 2.0D0*q(1) - 1.0D0
      dz(2) = 2.0D0*q(2) - 1.0D0

      end subroutine

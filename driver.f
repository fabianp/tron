      program tron

      integer nmax, nnzmax
      parameter(nmax=90000, nnzmax=10*nmax)
      integer nread, nwrite
      parameter(nread=1,nwrite=2)
c     **********
c
c     Driver for bound-constrained problems.
c
c     Subprograms called
c
c       USER ........... dminfg, dminhs, dminsp, dminxb
c
c       MINPACK-2 ...... dtron, dgpnrm2, dsphesd, dmid, wallclock
c
c       Level 1 BLAS ... dnrm2
c
c     MINPACK-2 Project. March 2000.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      double precision zero, one
      parameter(zero=0.0d0,one=1.0d0)

      character*60 task

      logical search
      integer arow_ind(nnzmax), acol_ind(nnzmax), brow_ind(nnzmax)
      integer arow_ptr(nmax+1), acol_ptr(nmax+1), bcol_ptr(nmax+1), 
     +        lrow_ind(nnzmax), lcol_ptr(nmax+1)
      integer indfree(nmax), isave(3), iwa(6*nmax) 
      integer n, itermax
      double precision f, delta
      double precision x(nmax), xl(nmax), xu(nmax), g(nmax) 
      double precision s(nmax), xc(nmax)
      double precision adiag(nmax), bdiag(nmax), ldiag(nmax)
      double precision a(nnzmax), b(nnzmax), l(nnzmax)
      double precision wa1(8*nmax), dsave(3)

c     Tolerances.

      double precision cgtol, frtol, fatol, fmin

c     Summary information.

      integer iterscg, nbind, nfev, nfree, ngev, nhev, nhsev

c     Evaluation of the Hessian matrix.

      integer maxgrp, numgrp
      integer listp(nmax), ngrp(nmax)
      double precision y(nmax), eta(nmax)

c     Test problems.

      character*6 prob
      integer i, info, j, maxfev, nnz, nx, ny
      double precision gnorm, gnorm0, par, gtol
      double precision dgpnrm2, dnrm2

c     Timing.

      double precision wctime1, wctime2, time1, time2
      double precision fgtime, htime, wctime

      open (nread,file='tron.dat',status='old') 
      open (nwrite,file='tron.info')

      do while (1. eq. 1)

c        Read problem data.
      
         read (nread,*) prob, n, nx, ny, par
         write (*,*)    prob, n, nx, ny, par
         
         if (prob(1:4) .eq. 'STOP') then
            close(nread)
            stop
         end if

c        Generate the starting point and the initial point.
         
         call dminfg(n,nx,ny,x,f,g,'XS',prob,par,wa1)
         call dminxb(n,nx,ny,xl,xu,prob)

c        Project the initial point into [xl,xu].

         call dmid(n,x,xl,xu)

c        Initialize variables.

         nfev = 0 
         ngev = 0
         nhev = 0
         nhsev = 0
         fgtime = zero
         htime = zero 

c        Set parameters. 
         
         itermax = n 
         maxfev = 1000
         fatol = zero
         frtol = 1.d-12
         fmin = -1.0d+32
         cgtol = 0.1d0
         gtol = 1.0d-5 
         
         call wallclock(wctime1)

c        Calculate the sparsity pattern.

         call dminsp(n,nx,ny,nnz,arow_ind,acol_ind,prob)

c        Transform coordinate format into compressed format.

         call dsetsp(n,nnz,arow_ind,acol_ind,acol_ptr,arow_ptr,1,
     +               info,listp,ngrp,maxgrp,iwa)
         if (info .le. 0) then
            write (nwrite,*) 'ERROR: INFO IN SUBROUTINE SETSP IS ',info
            stop
         end if

c        Information at starting point.
         
         nfree = 0
         nbind = 0
         do i = 1, n
            if (xl(i) .lt. x(i) .and. x(i) .lt. xu(i)) then
               nfree = nfree + 1
            else if ((x(i) .eq. xl(i) .and. g(i) .ge. zero) .or.
     +               (x(i) .eq. xu(i) .and. g(i) .le. zero) .or.
     +               (xl(i) .eq. xu(i))) then
               nbind = nbind + 1
            end if
         end do

c        Start the iteration.

         write (nwrite,1000) prob, n, maxgrp, dble(nnz)/n, nfree,
     +      n-nfree, nbind
         task = 'START'
         search = .true.

         do while (search)

c           Function evaluation.

            if (task .eq. 'F' .or. task .eq. 'START') then
               call wallclock(time1)
               call dminfg(n,nx,ny,x,f,g,'F',prob,par,wa1)
               call wallclock(time2)
               nfev = nfev + 1
               fgtime = fgtime + (time2 - time1)
            end if

c           Evaluate the gradient and the Hessian matrix.

            if (task .eq. 'GH' .or. task .eq. 'START') then
               call wallclock(time1)
               call dminfg(n,nx,ny,x,f,g,'G',prob,par,wa1)
               call wallclock(time2)
               ngev = ngev + 1
               fgtime = fgtime + (time2 - time1)

c              Evaluate the Hessian matrix.

               call wallclock(time1)
               do i = 1, n
                  s(i) = zero
                  eta(i) = one
               end do
               do numgrp = 1, maxgrp
                  do j = 1, n
                     if (ngrp(j) .eq. numgrp) s(j) = one
                  end do
                  call dminhs(n,nx,ny,x,s,y,prob,par,wa1)
                  call dsphesd(n,arow_ind,acol_ind,arow_ptr,acol_ptr,
     +                       listp,ngrp,maxgrp,numgrp,eta,y,a,adiag,iwa)
                  do j = 1, n
                     if (ngrp(j) .eq. numgrp) s(j) = zero
                  end do
               end do
               nhev = nhev + 1
               nhsev = nhsev + maxgrp
               call wallclock(time2)
               htime = htime + (time2 - time1)

            end if

c           Initialize the trust region bound.

            if (task .eq. 'START') then
               gnorm0 = dnrm2(n,g,1) 
               delta = dnrm2(n,g,1)
            end if

c           Call the optimizer.

            if (search) then
               call dtron(n,x,xl,xu,f,g,a,adiag,acol_ptr,arow_ind,
     +            frtol,fatol,fmin,cgtol,itermax,delta,task,
     +            b,bdiag,bcol_ptr,brow_ind,
     +            l,ldiag,lcol_ptr,lrow_ind,
     +            xc,s,indfree,
     +            isave,dsave,wa1,iwa)
            end if

c           Test for convergence or termination.
          
            if (task .eq. 'NEWX') then
               gnorm = dgpnrm2(n,x,xl,xu,g)
               if (gnorm .le. gtol*gnorm0) then
                  task = 'CONVERGENCE: GTOL TEST SATISFIED'
               end if
            end if
            if (nfev .gt. maxfev) then
               search = .false.
               task = 'ERROR: NFEV > MAXFEV'
            end if
            if (task(1:4) .eq. 'CONV') search = .false.

         end do

c        End of search.

         call wallclock(wctime2)
         write (*,*) task

c        Summary information.

         nfree = 0
         nbind = 0
         do i = 1, n
            if (xl(i) .lt. x(i) .and. x(i) .lt. xu(i)) then
               nfree = nfree + 1
            else if ((x(i) .eq. xl(i) .and. g(i) .ge. zero) .or.
     +               (x(i) .eq. xu(i) .and. g(i) .le. zero) .or.
     +               (xl(i) .eq. xu(i))) then
               nbind = nbind + 1
            end if
         end do
         
         iterscg =  isave(3)
         call dminfg(n,nx,ny,x,f,g,'FG',prob,par,wa1)
         gnorm = dgpnrm2(n,x,xl,xu,g)
         write (nwrite,2000) nfree, n-nfree, nbind, 
     +      nfev, ngev, nhev, nhsev, iterscg, f, gnorm

c        Timing information.
         
         wctime = wctime2 - wctime1
         fgtime = 100*(fgtime/wctime)
         htime  = 100*(htime/wctime)

         write (nwrite,3000) wctime, fgtime, htime, task

      end do

 1000 format (' Problem ', a6,                                   //,
     +        ' Number of variables                         ',i12,/,
     +        ' Number of coloring groups                   ',i12,/,
     +        ' Average number of nonzeros in the strictly  '    ,/,
     +        ' lower triangular part of the Hessian matrix ',f12.2,/,
     +        ' Number of free variables at x(start)        ',i12,/,
     +        ' Number of active variables at x(start)      ',i12,/,
     +        ' Number of binding variables at x(start)     ',i12,/)

 2000 format (
     +        ' Number of free variables at x(final)        ',i12,/,
     +        ' Number of active variables at x(final)      ',i12,/,
     +        ' Number of binding variables at x(final)     ',i12,/,
     +        ' Number of function evaluations              ',i12,/,
     +        ' Number of gradient evaluations              ',i12,/,
     +        ' Number of Hessian evaluations               ',i12,/,
     +        ' Number of Hessian-vector evaluations        ',i12,/,
     +        ' Number of conjugate gradient iterations     ',i12,//,
     +        ' Function value at final iterate          '   ,d15.8,/,
     +        ' Projected gradient at final iterate      '   ,d15.3,/)

 3000 format (' Total execution time                        ',f12.2,/,
     +        ' Percentage in function evaluations          ',f12.2,/,
     +        ' Percentage in Hessian evaluations           ',f12.2,//,
     +        ' Exit message     '                           ,a60,/)

      end


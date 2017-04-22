!! ECE697NA
!! Name: Poorna CHandra, Fathima
!!
!! This program proposes to read a matrix A in mtx (matrix-market/ coordinate), transform it
!! in CSR format and solve the linear Ax=b (elements are b are set equal to 1)
!! 
!! The code can be compiled as follows: "ifort -o hw5 hw5.f90 -mkl" using intel Fortran compiler and MKL
!! and example of runs: "./hw5 laplace3d 3 F"
!! the first argument is the name of the matrix file "laplace3d.mtx"
!! the second argument is the type of algo (see above)
!! the third argument is asking if all the entries are supplied (F)
!!  or only the lower (L) or upper part (U) for symmetric matrix 
!!
!! Remark: at the end of the program various utility routines are provided:
!! 
!! *dcsrmm: multiply csr matrix by vectors
!! *dcsrsv: solve lower/upper triangular system using csr matrix
!! *csrbwd: return bandwidth of a csr matrix
!! *dgbalu: perform LU banded factorization
!! *dgbalu2: needed for dgbalu
!! *dtbsm: solver lower/upper triangular system using banded matrix
!! *dcsr2csr_up: convert full csr to upper csr
!! *dcoo2csr: to convert coo to csr format

program spike_lu

!!!!!!! variable declaration
  implicit none
  integer :: i,j,N,k,Nx,Ny,e,nnz,rS,cS, Nj,M, s_cap_size
  integer :: t1,t2,tim
  double precision :: L,nres,err
  double precision,dimension(:),allocatable :: sa,b,r,ref_x,x,G,Gj, g_cap,x_cap
  integer,dimension(:),allocatable :: isa,jsa,ipiv
  character(len=100) :: name
  character(len=1) :: proc, uplo
  integer :: P
  !! for banded solver
  double precision,dimension(:,:),allocatable :: ba, ba_lu, Aj, spk_BjCj, spk_cap_mat, Bj, Cj
  double precision :: nzero,norm
  integer :: kl,ku,info
  !!for pardiso
  double precision,dimension(:),allocatable :: a
  integer,dimension(:),allocatable :: ia,ja
  character(len=1) :: cc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! read command line argument name="1","2","3","4" or "5"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,name)
  call getarg(2,proc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Read matrix "name" in matrix format                 !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  read(proc, '(i)') p
  print *, "Number of Processes:",p

  open(10,file=trim(name)//'.mtx',status='old')
  k=0
  cc='%'
  do while(cc=='%')
     k=k+1
     read(10,'(A1)') cc
  end do
  close(10)

  open(10,file=trim(name)//'.mtx',status='old')
  do i=1,k-1
     read(10,'(A1)') cc
  enddo
  read(10,*) n,n,nnz

  allocate(ia(nnz))
  allocate(ja(nnz))
  allocate(a(nnz))
  do i=1,nnz
     read(10,*) ia(i),ja(i),a(i)
  end do
  close(10)

!!!!!!!!!!! transform into csr format
  allocate(isa(1:N+1))
  allocate(jsa(1:nnz))
  allocate(sa(1:nnz))  
  call dcoo2csr(n,nnz,ia,ja,a,isa,jsa,sa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  print *,'Matrix:  ',trim(name)
  print *,'N=',N
  print *,'nnz=',nnz
  print *,''
!!!!!!!!!!!!!!


  allocate(b(1:N))! ! Right hand-side
  b=1.0d0

  allocate(ref_x(1:N)) !! solution
  ref_x=0.0d0 ! initial starting vector


  allocate(r(1:N)) !! residual array

     uplo = 'F'
     if (uplo/='F') then
        print *,'This version will not work if not in full sparse format'
        stop
     end if


!!! find the bandwidth
     call csrbwd(N,jsa,isa,kl,ku)
     print *,'banded format kl,ku',kl,ku

     if (kl+ku+1>1000) then
        print *,'bandwidth probably too large for banded storage'
        stop
     end if


!!! form banded matrix
     allocate(ba(1:kl+ku+1,1:N))
     ba=0.0d0 ! all elements at zero
     do i=1,N
        do k=isa(i),isa(i+1)-1
           ! (Row/Column) csr (i,jsa(k)) to banded (ku+1+i-jsa(k),jsa(k)
           ba(ku+1+i-jsa(k),jsa(k))=sa(k)
        enddo
     enddo

     !print *, "Reference Solution - The banded matrix is:", ba
!Reference Solution for Spike algorithms output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  Direct banded solver (SPIKE banded primitives)!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call system_clock(t1,tim) ! initialize time
!!! Solve Ax=b

     ! LU factorization
     nzero=0.0d0
     norm=0.0d0
     call DGBALU(N,kl,ku,ba,kl+ku+1,nzero,norm,info)
     !print *, "Reference Solution - The LU banded matrix is:", ba
     print *,'info',info
     ! LU solve (2 steps)
     ref_x=b
     call DTBSM('L','N','U',N,1,kl,bA(ku+1,1),kl+ku+1,ref_x,N)
     call DTBSM('U','N','N',N,1,ku,bA,kl+ku+1,ref_x,N)
     print *, "Reference Solution is:", ref_x

     ! compute residual
     r=b
     call dcsrmm(uplo,N,N,1,-1.0d0,sa,isa,jsa,ref_x,1.0d0,r)
     nres=sum(abs(r))/sum(abs(b)) ! norm relative residual 
     print *, "Residual from nres is", nres


  call system_clock(t2,tim) ! final time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *,'Total time',(t2-t1)*1.0d0/tim
! Refernce Solution calculation ends

!! Recreating banded A matrix again - while generating reference solution ba was modified by DGBALU function!!
     ba=0.0d0 ! all elements at zero
     do i=1,N
        do k=isa(i),isa(i+1)-1
           ! (Row/Column) csr (i,jsa(k)) to banded (ku+1+i-jsa(k),jsa(k)
           ba(ku+1+i-jsa(k),jsa(k))=sa(k)
        enddo
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pre-processing stage:: partitioning phase!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Nj = N/P
M = (kl+ku+1-1)/2
print *, "Order of each diagonal blocks - Aj: ", Nj  !! This needs to be put in variable still
print *, "Order of each sub diagonal blocks - Bj, Cj: ", M  !! This needs to be put in variable still
!and also this can change if bandwidth ku ad kl are equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pre-processing stage:: Factorization phase!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!Extracting Block Diagonal Matrix
!! No need to explicitly extract Diagonal Elements
!! It will be part of Matrix ba_lu - which is computed below; each Aj in it is in LU format

!!Generating Reduced Spike Matrix
! No need to generate Spike Matrix S as it is not needed
 
 s_cap_size = 2*(P-1)*M;
 
 nzero=0.0d0
 norm=0.0d0
 allocate(Aj(kl+ku+1,Nj))
 allocate(spk_BjCj(Nj,M))
 allocate(Bj(M,M))
 allocate(Cj(M,M))
 allocate(spk_cap_mat(s_cap_size, s_cap_size))
 allocate(ba_lu(1:kl+ku+1,1:N))
 
 spk_cap_mat = 0.0d0
 ba_lu = 0.0d0
 spk_BjCj = 0.0d0
 Bj = 0.0d0
 Cj = 0.0d0

 do k=1, P
 !!LU Factorization of A
   Aj = ba(1:kl+ku+1,(1+(k-1)*Nj):(k*Nj))
   nzero=0.0d0
   norm=0.0d0
   print *, "K is",k
   call DGBALU(Nj,kl,ku,Aj,kl+ku+1,nzero,norm,info)

   print *,'info',info
   ba_lu(1:kl+ku+1,(1+(k-1)*Nj):(k*Nj)) = Aj

   !!Fetching Bj from next partition - spk_BjCj is acting as spk_Bj here
   if(k<P) then
     spk_BjCj = 0.0d0
     Bj = 0.0d0
     !call DLACPY('L', ku, ku, ba(1,(1+(k)*Nj)), kl+ku, Bj, ku) !! When in doubt, check this :p
     call DLACPY('L', ku, ku, ba(1,(1+(k)*Nj)), kl+ku, spk_BjCj(Nj-M+1:Nj, 1:M), ku) !! When in doubt, check this :p
    ! print *, "Spk Bj is", spk_BjCj
     !print *, "Bj is", Bj

     !LU SOlve
     call DTBSM('L','N','U',Nj,M,kl,Aj(ku+1,1),kl+ku+1,spk_BjCj,Nj)
     call DTBSM('U','N','N',Nj,M,ku,Aj,kl+ku+1,spk_BjCj,Nj)

     !Storing Vj bottom in reduced spike matrix - spk_BjCj is Vj now
     spk_cap_mat(1+(k-1)*2*M : 1+(k-1)*2*M+(M-1), (2*k*M)-(M-1) : 2*k*M) = spk_BjCj(Nj-M+1:Nj, 1:M)
     !Storing Vj top in reduced spike matrix - spk_BjCj is Vj now
     if(k>1) then
       spk_cap_mat((2*k-3)*M+1 : (2*k-2)*M, (2*k*M)-(M-1) : 2*k*M) = spk_BjCj(1:M, 1:M)
     end if
   end if

   !!Fetching Cj from previous partition - spk_BjCj is acting as spk_Cj here
   if(k>1) then
     spk_BjCj = 0.0d0
     Aj = ba_lu(1:kl+ku+1,(1+(k-1)*Nj):(k*Nj)) !Current partitiion Aj
     !print *, "Aj is", Aj
     !call DLACPY('U', kl, kl, ba(ku+1 : kl+ku+1,(1+(k-2)*Nj+Nj-kl):((k-1)*Nj)), kl+ku, Cj, kl) !! When in doubt, check this :p
     !call DLACPY('U', kl, kl, ba(kl+ku+1,(1+(k-2)*Nj+Nj-kl)), kl+ku, Cj, kl) !! When in doubt, check this :p
     call DLACPY('U', kl, kl, ba(kl+ku+1,(1+(k-2)*Nj+Nj-kl)), kl+ku, spk_BjCj(1:M, 1:M), ku) !! When in doubt, check this :p
     !print *, "Cj is", spk_BjCj

     !LU SOlve - Wj is nothing spk_BjCj here
     call DTBSM('L','N','U',Nj,M,kl,Aj(ku+1,1),kl+ku+1,spk_BjCj,Nj)
     call DTBSM('U','N','N',Nj,M,ku,Aj,kl+ku+1,spk_BjCj,Nj)

     !Storing Wj top in reduced spike matrix - spk_BjCj is Wj now
     spk_cap_mat((2*(k-1)*m)-(M-1) : 2*(k-1)*M, 1+(k-2)*2*M : 1+(k-2)*2*M+(M-1)) = spk_BjCj(1:M, 1:M)
     !Storing Wj bottom in reduced spike matrix - spk_BjCj is Wj now
     if(k<P) then
       spk_cap_mat((2*(k-1)*M)+1 : (2*k-1)*M, 1+(k-2)*2*M : 1+(k-2)*2*M+(M-1)) = spk_BjCj(Nj-M+1:Nj, 1:M)
     end if
   end if
 enddo

 do i = 1, s_cap_size
   spk_cap_mat(i,i) = 1.0d0
 enddo
 !! Spike Matrix is ready!!. Hopefully :p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Post-Processing stage:: D * G = F        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This serialy solving of G needs to be parallelized later

!Note: Aj is used here instead of "D" Diagonal block matrix assuming 
!DGBALU refers to only the diagonal block

 allocate(G(N))
 G = 0.0d0
 allocate(Gj(Nj))
 Gj = 0.0d0
 do k=1,P
   !!In the below step, using already existing LU factorized "A" Matrix
   Aj = ba_lu(1:kl+ku+1,(1+(k-1)*Nj):(k*Nj))
   !LU SOlve
   Gj = b(1+(k-1)*Nj:k*Nj)
   !print *, "Gj is", Gj
   call DTBSM('L','N','U',Nj,1,kl,Aj(ku+1,1),kl+ku+1,Gj,Nj)
   !print *, "Gj is", Gj
   call DTBSM('U','N','N',Nj,1,ku,Aj,kl+ku+1,Gj,Nj)
   !print *, "Gj is", Gj
   G(1+(k-1)*Nj:k*Nj) = Gj

 enddo
 !print *, "Final G is", G

 !! Form G_cap
 allocate(g_cap(s_cap_size))

 do i = 1,P
   !To store Gb
   if(i<P) then
     g_cap(1+(i-1)*2*M : (i-1)*2*M+M) = G((i*Nj)-M+1 : (i*Nj))
   end if
   !To store Gt
   if(i>1) then
     g_cap(m*((2*i)-3)+1 : 2*M*(i-1)) = G(((i-1)*Nj)+1 : ((i-1)*Nj+M))
   end if
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Post-Processing stage:: S * X = G        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This is just for reference; later this needs to be replaced with solution for
!! reduced form of S, X and G

  allocate(x(1:N)) !! solution
  x=0.0d0 ! initial starting vector
  allocate(x_cap(1:s_cap_size)) !! solution
  x_cap=0.0d0 ! initial starting vector
  allocate(ipiv(N))
  call DGETRF(s_cap_size,s_cap_size,spk_cap_mat,s_cap_size,ipiv,info)
  x_cap=g_cap
  !print *, "G is", G

  call DGETRS('N', s_cap_size, 1, spk_cap_mat, s_cap_size, ipiv, x_cap, s_cap_size, info)
  print *,'info',info

  print *, "Final Solution is:",x_cap

  !r=b
  !call DGEMM('N', 'N', N, 1, N, 1.0d0, A_mat, N, x, N, -1.0d0, r, N)
  !print *, "Residual from maxval is", maxval(abs(r))/maxval(abs(b))
  !nres=sum(abs(r))/sum(abs(b)) ! norm relative residual 
  !print *, "Residual from nres is", nres

end program spike_lu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dcsrmm(UPLO,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*A*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  implicit none
  character(len=1) :: UPLO
  integer :: N,M,rhs
  double precision :: alpha,beta
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(M,*):: x
  double precision,dimension(N,*) ::b
!!!!!!!!        

  integer ::i,j,k,s
  double precision,parameter :: DZERO=0.0d0
  double precision,dimension(rhs) :: summ

  if (UPLO=='F') then

     do i=1,N
        summ=DZERO
        do k=ia(i),ia(i+1)-1
           j=ja(k)
           summ(1:rhs)=summ(1:rhs)+alpha*a(k)*x(j,1:rhs)
        end do
        b(i,1:rhs)=beta*b(i,1:rhs)+summ(1:rhs)
     end do

  elseif ((UPLO=='L').or.(UPLO=='U')) then

     do s=1,rhs

        if (beta/=DZERO) then
           b(1:N,s)=beta*b(1:N,s)
        else
           b(1:N,s)=DZERO
        endif


        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,s)=b(i,s)+alpha*a(k)*x(j,s)
              if (j/=i) b(j,s)=b(j,s)+alpha*a(k)*x(i,s)               
           end do
        end do

     enddo

  end if

end SUBROUTINE dcsrmm





SUBROUTINE dcsrsv(UPLO,N,rhs,a,ia,ja,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine solve
  !
  !      b=A^-1*x      (Ax=b A=L or A=U)
  !  
  !      where A is a NxN triangular matrix stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='L'  A is lower triangular
  !       if UPLO='U'  A is upper triangular 
  !
  !  N    (input) INTEGER
  !        The number of row/column of the matrix A and row of matrix B.  N >= 0.
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  implicit none
  character(len=1) :: UPLO
  integer :: N,rhs
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(N,*) ::b
!!!!!!!!        

  integer ::i,j,k
  double precision,dimension(rhs) :: summ



  if (UPLO=='L') then


     DO  i=1,N
        summ=b(i,1:rhs)
        DO  k=ia(i), ia(i+1)-1
           if (ja(k)<i)  summ=summ-a(k)*b(ja(k),1:rhs)
           if (ja(k)==i) j=k
        end DO
        b(i,1:rhs)=summ/a(j)
     end DO

  elseif (UPLO=='U') then

     DO  i=N,1,-1
        summ=b(i,1:rhs)
        DO  k=ia(i), ia(i+1)-1
           if (ja(k)>i)  summ=summ-a(k)*b(ja(k),1:rhs)
           if (ja(k)==i) j=k
        end DO
        b(i,1:rhs)=summ/a(j)
     end DO



  end if


end SUBROUTINE dcsrsv





subroutine csrbwd(n,ja,ia,kl,ku)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Purpose: computes the lower and upper bandwidth of a sparse(CSR) matrix 
  ! 
  ! Arguments: 
  !       n (input) integer, size of the matrix 
  !       ja (input) integer, column pointers 
  !       ia (input) integer, row indices 
  !       kl (output) intger, lower bandwidth 
  !       ku (output) integer, upper bandwidth 

  integer,dimension(*) :: ja
  integer  :: n
  integer,dimension(1:n+1) ::ia
!!!!!!!!!!!!
  integer ::  kl,ku
  integer :: dst,i,k

  kl = - n
  ku = - n
  do  i=1,n
     do  k=ia(i),ia(i+1)-1
        dst = i-ja(k)
        kl = max(kl,dst)
        ku = max(ku,-dst)
     end do
  end do
end subroutine csrbwd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine DGBALU(N,kl,ku,A,LDA,nzero,norm,info)
  !  Purpose
  !  =======
  !
  !  it computes an ("approximate") LU factorization (blocked-BLAS3) of a real   n-by-n band matrix A
  !  without using partial pivoting with row interchanges.+ diagonal boosting if any
  !
  !  This is the blocked version of the algorithm for square matrices, calling Level 3 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  N      (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !  KU     (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  A      (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
  !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  !
  !          On exit, contains the LU factorization: U is stored as an
  !          upper triangular band matrix with KU superdiagonals in
  !          rows 1 to KU+1. L is stored as an lower triangular band matrix
  !          with KL subdiagonals in rows KU+2 to KL+KU+1. The unit diagonal of L
  !          is not stored.
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
  !           a boost will take place if pivot <nzero*norm
  !
  ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
  !          Value of the boosting will be (+/-) eps*norm
  !          norm is a dummy argument (unused) if nzero=0.0d0
  ! INFO    (output) INTEGER
  !          = 0: successful exit
  !          -6 < INFO < 0 : if INFO = -i, the i-th argument had an illegal value
  !          =-10: internal memory allocation problem
  !          > 0: if INFO = +i, The factorization has been completed and 
  !               +i "nzero" pivots have been found (so +i boost took place) 
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================

  implicit none
  integer :: N,kl,ku,LDA
  double precision,dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info

!!!!!!!!!!!!!!!!!!
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  integer :: m,i,nbm,jmin,klr,kur,j,k,il,ju,shift,info2
  double precision,dimension(:,:),pointer :: taux1,taux2
  INTEGER:: ILAENV
  integer :: info_alloc
!!!!!!!!!!!!!!!!!!!
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( N.LT.0 ) THEN
     INFO = -1
  ELSE IF( KL.LT.0 ) THEN
     INFO = -2
  ELSE IF( KU.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGBALU', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN

  !     Determine the block size for this environment
  m = ILAENV( 1, 'DGBTRF', ' ', N, N, KL, KU )
  m = MIN( m, min(kl,ku) )

  if (m==1) then
     !! use unblocked algorithm

     CALL DGBALU2(n,n,KL,KU,A,LDA,nzero,norm,info2)
     if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'DGBALU', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

  else

     !! use block algorithm
     nbm=n/m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Start recursion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! we consider the (sub)matrix
!!! A B
!!! C D
     ! we want (M1 is also the end of L1---- LU is then done on a rectangular matrix)
!!! L1       U1 K1 
!!! M1 L2       U2    
!!!
!!!!!!!! M1 (resp. K1) may be composed by rectangular part and triangular parts.

     allocate(taux1(1:m, 1:m)) 
     allocate(taux2(1:m, 1:m)) 


     taux1=DZERO
     taux2=DZERO

     if (n-1<=max(kl,ku)) jmin=-m+1

     DO  k=1,nbm-1

        if (n-(k)*m-1>max(kl,ku)) then

           jmin=(k-1)*m+1

           !For M1 (rectangular remaining row)
           klr=min(kl,n-jmin)-m
           if (klr<0) klr=0
           !For K1 (rectangular remaining column)
           kur=min(ku,n-jmin)-m
           if (kur<0) kur=0
           !! total row for submatrix -- L1 (including M1)
           il=min(kl,n-(jmin+m-1))+m

!!!!!!! L1U1 factorization of the rectangular (il*m) block A

           CALL  DGBALU2(il,m,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
           if (info2<0) THEN
              info=info2+1
              if (INFO.NE.0) CALL XERBLA( 'DGBALU', -INFO )
              RETURN
           ELSE
              info=info+info2
           END IF

!!!!!
!!!!!!!!  L1K1=B so find K1
!!!!!!

           !! B may be composed by a rectangular part (m*kur) and triangular part (m*m)  
           ! rectangular part (update A directly)
           if (kur>0) CALL DTRSM('L','L','N','U',m,kur,DONE,A(ku+1,jmin),kl+ku,A(ku+1-m,jmin+m),kl+ku) 
           ! triangular part (update A indirectly)   
           do i=1,m
              do j=1,i
                 taux1(i,j)=A(ku+i-m-kur+1-j,jmin+m+kur-1+j) 
              end do
              do j=i+1,m
                 taux1(i,j)=DZERO
              end do
           enddo

           CALL DTRSM('L','L','N','U',m,m,DONE,A(ku+1,jmin),kl+ku,taux1,m) 
           do j=1,m
              do i=j,m
                 A(ku+1+i-j-m-kur,jmin+m+kur+(j-1))=taux1(i,j) 
              end do
           end do

!!!!!!
!!!!!!!!  UPDATE D to close the recursion D<=D-M1K1
!!!!!!
           do i=1,m
              do j=i,m
                 taux2(i,j)=A(ku+m+klr+i+1-j,jmin-1+j)
              enddo
           enddo


!!!! we can decompose here into 4 multiplication

           ! M1 (rectangular klr*m) * K1 (rectangular m*kur)
           if ((klr>0).and.(kur>0)) call DGEMM('N','N',klr,kur,m,-DONE,A(ku+1+m,jmin),kl+ku,A(ku+1-m,jmin+m),&
                &kl+ku,DONE,A(ku+1,jmin+m),kl+ku)
           ! M1 (rectangular klr*m) * K1 (lower triangular m*m)
           if (klr>0) then
              call DGEMM('N','N',klr,m,m,-DONE,A(ku+1+m,jmin),kl+ku,taux1,m,DONE,A(ku+1-kur,jmin+m+kur),kl+ku)
           end if
           ! M1 (upper triangular m*m) * K1 (rectangular m*kur)
           if (kur>0) then
              call DGEMM('N','N',m,kur,m,-DONE,taux2,m,A(ku+1-m,jmin+m),kl+ku,DONE,A(ku+1+klr,jmin+m),kl+ku)
           end if
           ! M1 (upper triangular m*m) * K1 (lower triangular m*m)
           if (klr>0) then
              if (kur>0) then 
                 shift=-kur+klr
              else
                 shift=klr
              endif
           else
              if (kur>0) then
                 shift=-kur
              else
                 shift=0
              endif
           endif

           call DGEMM('N','N',m,m,m,-DONE,taux2,m,taux1,m,DONE,A(ku+1+shift,jmin+m+kur),kl+ku)
        end if
     end DO

!!!!!!!!!! LU on the Last block

     jmin=jmin+m
     il=n-jmin+1
     ju=il

!!!!!!! L1U1 factorization of the rectangular (il*ju) block A
     CALL  DGBALU2(il,ju,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
     if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'DGBALU', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

     deallocate(taux1)  
     deallocate(taux2)  

  end if

end subroutine DGBALU






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DGBALU2(M,N,kl,ku,A,LDA,nzero,norm,info)
  !
  !  Purpose
  !  =======
  !
  !  it computes an ("approximate") LU factorization (unblocked-BLAS2) of a real  m-by-n band matrix A
  !  without using partial pivoting with row interchanges.+ diagonal boosting if any
  !
  !  This is the unblocked version of the algorithm, calling Level 2 BLAS, inspired from LAPACK DGBTF2.
  !
  !  Arguments
  !  =========
  !
  !  M      (input) INTEGER
  !          The number of rows of the matrix A.  M >= 0.
  !  N      (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !  KU     (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  A      (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
  !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
  !           a boost will take place if pivot <nzero*norm
  !
  ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
  !          Value of the boosting will be (+/-) eps*norm
  !          norm is a dummy argument (unused) if nzero=0.0
  ! INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -i, the i-th argument had an illegal value
  !          > 0: if INFO = +i, The factorization has been completed and 
  !               +i "nzero" pivots have been found (so +i boost took place) 
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  integer :: M,N,kl,ku,LDA
  double precision,dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info

!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  double precision :: pivot
  integer :: j,JU,KM
!!!!!!!!!!!!!!!!

  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSEIF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGBALU2', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN
  IF (nzero.EQ.DZERO) then
     pivot=DZERO
  else
     pivot=nzero*norm
  end IF
  JU=1
  DO  J = 1,MIN(M,N)
     if (abs(A(KU+1,J))<=pivot) then
        if (nzero.EQ.DZERO) then 
           info=-7
           CALL XERBLA( 'DGBALU2', -INFO )
           RETURN
        end IF
        A(KU+1, J)=A(KU+1,J)+sign(nzero,A(KU+1,J))*norm
        info=info+1
     end if
     KM=min(KL,M-J)
     JU=max(JU,min(J+KU,N))
     IF( KM.GT.0 ) THEN
        CALL DSCAL(KM, DONE/A(KU+1, J), A(KU+2, J), 1)
        IF( JU.GT.J )   CALL DGER(KM, JU-J, -DONE, A(KU+2, J), 1,A(KU, J+1), KL+KU, A(KU+1, J+1), KL+KU)
     end IF
  END DO

end subroutine DGBALU2







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DTBSM(UPLO,TRANS,DIAG,N,rhs,kd,A,LDA,B,LDB)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Purpose:solves one of the systems of equations in real double precision
  !
  !     A*x = b,   or   A'*x = b,  
  !
  !  where b and x are n by rhs element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular band matrix, with ( kd + 1)
  !  diagonals.
  !
  !  It is a  Blocked version (BLAS3) of the routine DTBSV which allow blocking with 
  !  multiple right-hand sides      
  !
  ! Arguments: 
  !       UPLO (input) character: if 'L' or 'l' lower triangular , if 'U' or 'u' upper triangular
  !       TRANS (input) character:specifies the equations to be solved as follows:
  !              TRANS = 'N' or 'n'   A*x = b.
  !              TRANS = 'T' or 't'   A'*x = b.
  !       DIAG (input) character: if 'N' or 'n' non-unit diagonal, if 'U' or 'u' unit diagonal
  !       N (input) integer: order of matrix A
  !       rhs (input) integer :number of right-hand-sides for vector B        
  !       kd (input) integer: lower bandwidth of the matrix A if UPLO='L', or
  !                           upper bandwidth of the matrix A if UPLO='U'
  !       A (input) double precision array, banded input matrix 
  !           Note: same format than the one used in BLAS DTBSV
  !
  !           Before entry with UPLO = 'U' or 'u', the leading ( kd + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( kd + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row kd, and so on. The top left kd by kd triangle
  !           of the array A is not referenced.
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( kd + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right kd by kd triangle of the
  !           array A is not referenced.
  !           
  !       LDA (input) integer : The leading dimension of the array A.  LDA >= kd+1.
  !       B (input/output) double precision 2D array: right hand side contains the solution on exit 
  !       LDB (input) integer : The leading dimension of the array B.  LDB >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  !======================================================================
  implicit none
  integer :: N,LDA,LDB,rhs
  DOUBLE PRECISION,dimension(LDA,*),intent(in) ::   A
  character(len=1)  :: UPLO,DIAG,TRANS
  integer :: kd
  DOUBLE PRECISION,dimension(LDB,*) ::   B

!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,parameter :: DONE=1.0d0
  character(len=1) :: TRANSB,SIDE,UPLO2,DIAG2
  integer :: k,ijinit,i
  integer :: m,m1,max
  integer :: i1,i2,i3,i4,i5,i6,j1,j2,j3,j4,s1,s2
  double precision,dimension(kd,rhs) :: aux
  integer :: ldaux
!!!!!!!!!!!!!!!!!!!!!!!!



  !! Size of the Recursive diagonal block and initialization of UPLO2
  !! for triangular off-diagonal blocks
  m=min(kd,n)


!!!!!!!!!!!!!! Case n<kd
  if (n<kd) then ! treat as dense
     if ((UPLO=='L').or.(UPLO=='l')) then
        call DTRSM('L',UPLO,TRANS,DIAG,N,rhs,DONE,A(1,1),LDA-1,B,LDB)
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        call DTRSM('L',UPLO,TRANS,DIAG,N,rhs,DONE,A(kd+1,1),LDA-1,B,LDB)
     endif
     return
  endif
!!!!!!!!!!!!


  if ((UPLO=='L').or.(UPLO=='l')) then
     UPLO2='U'
  elseif ((UPLO=='U').or.(UPLO=='u')) then
     UPLO2='L'
  end if

  !!initialization of DIAG2,SIDE
  DIAG2='N'
  SIDE='L'


  !! how many submatrices (m*m) in A
  max=n/kd
  !! size of the first diagonal block before recursion
  m1=n-max*kd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Initialisation before the loop
!!!! only once for the first block decomposition with size m1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (m1/=0) then !! if there is a first block

     ldaux=kd

     TRANSB = 'N'

     if ((TRANS=='N').or.(TRANS=='n')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=1
           i1=1
           j3=j1
           i3=i1+m1
           j4=j1
           i4=kd+1

           i5=j1+m1 !rhs
           i6=i5+(kd-m1)

        elseif((UPLO=='U').or.(UPLO=='u')) then
           j1=n-m1+1                
           i1=kd+1
           j3=j1
           i3=m1+1
           j4=j1
           i4=1

           i5=j1-(kd-m1)
           i6=i5-m1
        end if

     elseif ((TRANS=='T').or.(TRANS=='T')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=n-m1+1
           i1=1
           j3=n-kd+1
           i3=kd+1-m1  
           j4=j1-kd
           i4=kd+1

           i5=j1-(kd-m1)
           i6=i5-m1
        elseif ((UPLO=='U').or.(UPLO=='u')) then
           j1=1
           i1=kd+1
           j3=j1+m1
           i3=i1-m1
           j4=kd+1
           i4=1

           i5=j1+m1
           i6=i5+(kd-m1)
        end if
     end if

!!!!!! Solve using blas3 the initial block
     call DTRSM ( SIDE, UPLO, TRANS, DIAG, m1, rhs, DONE, A(i1,j1), LDA-1,B(j1,1), LDB )
!!!!!! update the right hand side using blas3 in 2 steps
     ! step 1
     call DGEMM ( TRANS, TRANSB, kd-m1, rhs, m1, -DONE, A(i3,j3), LDA-1, B(j1,1), LDB,DONE, B(i5,1), LDB ) 
     ! step 2
     do i=1,rhs
        call DCOPY(m1,B(j1,i),1,aux(1,i),1)
     end do
     call DTRMM ( SIDE, UPLO2, TRANS, DIAG2, m1, rhs, DONE, A(i4,j4), LDA-1,aux,ldaux ) 
     do i=1,rhs
        call DAXPY(m1,-DONE,aux(1,i),1,B(i6,i),1)
     end do

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! starting point origin !! new origin and other initialization before loop
  if ((TRANS=='N').or.(TRANS=='n')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=m1+1
        i1=1
        i2=kd+1
        s1=1
        s2=0
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=n-m1-kd+1
        i1=kd+1
        i2=1
        s1=-1
        s2=0
     end if
  elseif ((TRANS=='T').or.(TRANS=='t')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=n-m1-kd+1
        i1=1
        i2=kd+1
        s1=-1
        s2=-1
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=m1+1
        i1=kd+1
        i2=1
        s1=1
        s2=1
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LOOP for all the submatrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO k=1,max
     j1=ijinit+s1*(k-1)*m
     j2=j1+s2*kd
     i3=j1+s1*kd
     !! Solve using blas3 the initial bloc
     call DTRSM ( SIDE, UPLO, TRANS, DIAG, kd, rhs, DONE, A(i1,j1),LDA-1,B(j1,1), LDB ) 
     !! update the right hand side using blas3
     if (k/=max) then
        do i=1,rhs
           call DCOPY(kd,B(j1,i),1,aux(1,i),1)
        end do
        call DTRMM ( SIDE, UPLO2, TRANS, DIAG2, kd, rhs, DONE, A(i2,j2), LDA-1,aux,kd ) 
        do i=1,rhs
           call DAXPY(kd,-DONE,aux(1,i),1,B(i3,i),1)
        end do
     end if
  END DO



end subroutine DTBSM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dcsr2csr_up(opt,N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine extract the upper triangular part (including diagonal) of a real CSR matrix A 
  !      into a new CSR matrix sA
  !      the size of sa, and sja for csr upper format may be overestimated 
  !
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs extraction, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine returns only values of the array sia of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value sia(n+1) - 1 is the actual 
  !         number of the elements in the arrays sa and sja.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array sa and sja
  !
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k
  logical :: row,colval
!!!!
  row=.false.
  colval=.false.
  if (opt==0) then
     row=.true.
     colval=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     colval=.true.
  end if



  k=0
  if (row) sia(1)=1
  do i=1,n 
     if (row) sia(i+1)=sia(i)
     do j=ia(i),ia(i+1)-1
        if (ja(j)>=i) then
           k=k+1
           if (row) sia(i+1)=sia(i+1)+1
           if (colval) then
              sja(k)=ja(j)
              sa(k)=a(j)
           end if
        endif
     end do
  end do
end subroutine dcsr2csr_up






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  double precision,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  double precision :: adum


  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do



  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do


  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine dcoo2csr





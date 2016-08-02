module cons_project

contains

subroutine project(rho, u, v, w, rhoe, &
                   idir, iface, Scons, Scons_proj, version, notrace_in)

  ! perform a characteristic projection with tracing on vector S,
  ! resulting in vector Sproj here idir is the coordinate direction
  ! and iface = 0 for the lower face and 1 for the upper face.

  use meth_params_module, only: gamma_const, gamma_minus_1

  implicit none

  integer, parameter :: ncons = 5

  double precision , intent(in)  :: rho, u, v, w, rhoe
  integer          , intent(in)  :: idir, iface
  double precision , intent(in)  :: Scons(ncons)
  double precision , intent(out) :: Scons_proj(ncons)
  integer          , intent(in)  :: version
  logical, optional, intent(in)  :: notrace_in

  double precision :: l(ncons,ncons), r(ncons,ncons), e(ncons)
  double precision :: p, c
  double precision :: Vsq, H, A, g1

  integer :: n
  logical :: notrace

  if (present(notrace_in)) then
     notrace = notrace_in
  else
     notrace = .false.
  endif

  g1 = gamma_minus_1

  ! compute some primitive variables
  ! p = (rhoE - 0.5d0*rho*(u**2 + v**2 + w**2))*g1
  p = g1 * rhoe

  c = sqrt(gamma_const*p/rho)

  ! compute the eigenvalues
  select case (idir)
  case (1)
     e(1) = u - c
     e(2) = u
     e(3) = u
     e(4) = u
     e(5) = u + c

  case (2)
     e(1) = v - c
     e(2) = v
     e(3) = v
     e(4) = v
     e(5) = v + c

  case (3)
     e(1) = w - c
     e(2) = w
     e(3) = w
     e(4) = w
     e(5) = w + c
  end select

  ! compute the eigenvectors
  Vsq = u*u + v*v + w*w
  H = 0.5*Vsq + c**2/g1

  select case (idir)
  case (1)
     r(1,:) = [1.0d0, u - c, v,     w,     H - u*c  ]
     r(2,:) = [1.0d0, u,     v,     w,     0.5d0*Vsq]
     r(3,:) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, v        ]
     r(4,:) = [0.0d0, 0.0d0, 0.0d0, 1.0d0, w        ]
     r(5,:) = [1.0d0, u + c, v,     w,     H + u*c  ]

     A = 0.5*g1/c**2
     l(1,:) = A * [H + (c/g1)*(u-c),   -(u + c/g1), -v,            -w,            1.0d0 ]
     l(2,:) = A * [-2*H + (4/g1)*c**2, 2.0d0*u,     2.0d0*v,       2.0d0*w,       -2.0d0]
     l(3,:) = A * [-2.0d0*v*c**2/g1,   0.0d0,       2.0d0*c**2/g1, 0.0d0,         0.0d0 ]
     l(4,:) = A * [-2.0d0*w*c**2/g1,   0.0d0,       0.0d0,         2.0d0*c**2/g1, 0.0d0 ]
     l(5,:) = A * [H - (c/g1)*(u+c),   -u + c/g1,   -v,            -w,            1.0d0 ]

  case (2)
     r(1,:) = [1.0d0, u,     v - c, w,     H - v*c  ]
     r(2,:) = [0.0d0, 1.0d0, 0.0d0, 0.0d0, u        ]
     r(3,:) = [1.0d0, u,     v,     w,     0.5d0*Vsq]
     r(4,:) = [0.0d0, 0.0d0, 0.0d0, 1.0d0, w        ]
     r(5,:) = [1.0d0, u,     v + c, w,     H + v*c  ]

     A = 0.5*g1/c**2
     l(1,:) = A * [H + (c/g1)*(v-c),   -u,            -(v + c/g1), -w,            1.0d0 ]
     l(2,:) = A * [-2.0d0*u*c**2/g1,   2.0d0*c**2/g1, 0.0d0,       0.0d0,         0.0d0 ]
     l(3,:) = A * [-2*H + (4/g1)*c**2, 2.0d0*u,       2.0d0*v,     2.0d0*w,       -2.0d0]
     l(4,:) = A * [-2.0d0*w*c**2/g1,   0.0d0,         0.0d0,       2.0d0*c**2/g1, 0.0d0 ]
     l(5,:) = A * [H - (c/g1)*(v+c),   -u,            -v + c/g1,   -w,            1.0d0 ]

  case (3)
     r(1,:) = [1.0d0, u,     v,     w - c, H - w*c  ]
     r(2,:) = [0.0d0, 1.0d0, 0.0d0, 0.0d0, u        ]
     r(3,:) = [0.0d0, 0.0d0, 1.0d0, 0.0d0, v        ]
     r(4,:) = [1.0d0, u,     v,     w,     0.5d0*Vsq]
     r(5,:) = [1.0d0, u,     v,     w + c, H + w*c  ]

     A = 0.5*g1/c**2
     l(1,:) = A * [H + (c/g1)*(w-c),   -u ,           -v,            -(w + c/g1), 1.0d0 ]
     l(2,:) = A * [-2.0d0*u*c**2/g1,   2.0d0*c**2/g1, 0.0d0,         0.0d0,       0.0d0 ]
     l(3,:) = A * [-2.0d0*v*c**2/g1,   0.0d0,         2.0d0*c**2/g1, 0.0d0,       0.0d0 ]
     l(4,:) = A * [-2*H + (4/g1)*c**2, 2.0d0*u,       2.0d0*v,       2.0d0*w,     -2.0d0]
     l(5,:) = A * [H - (c/g1)*(w+c),   -u,            -v,            -w + c/g1,   1.0d0 ]

  end select


  ! Just copy Scons into Scons_proj
  if (notrace) then
     Scons_proj(:) = Scons(:)

  ! Project along characteristics
  else

     Scons_proj(:) = 0.0d0

     if (version .eq. 1) then

        if (iface == 0) then
           ! left moving
           do n = 1, ncons
              if (e(n) <= 0.0d0) then
                 Scons_proj(:) = Scons_proj(:) + dot_product(l(n,:), Scons(:))*r(n,:)
              endif
           enddo
   
        else
           ! right moving
           do n = 1, ncons
              if (e(n) >= 0.0d0) then
                 Scons_proj(:) = Scons_proj(:) + dot_product(l(n,:), Scons(:))*r(n,:)
              endif
           enddo
        endif

     else if (version .eq. 2) then

        ! Lo side
        if (iface == 0) then
           ! If any of the waves reach this face then set to full state
           ! e(1) = u-c
           if (e(1) .lt. 0.d0) then
                 Scons_proj(:) = Scons(:)
           end if

        ! Hi side
        else
           ! If any of the waves reach this face then set to full state
           ! e(5) = u+c
           if (e(5) .gt. 0.d0) then
                 Scons_proj(:) = Scons(:)
           end if
        endif
     else 
        call bl_abort('BAD VERSION IN PROJECT_S')
     endif
  endif


end subroutine project
end module cons_project

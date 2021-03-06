module fault_module

    implicit none
    real(kind=8), parameter :: LAT2METER = 111133.84012073893 ! from clawpack.geoclaw.data
    integer :: nsubfaults
    type subfault
      real(kind=8) :: width, depth, slip, longitude, rupture_time, rise_time
      real(kind=8) :: xcb(2)
    end type subfault
    type(subfault), allocatable :: subfaults(:)
    real(kind=8) :: center(2), theta, xcb(2), final_rupture_rise_time

contains

    subroutine load_fault(fname)

        implicit none

        character*12, intent(in) :: fname

        integer :: i
        real(kind=8) :: input_line(12), xp1, yp1, xp2, yp2, total_width

        call opendatafile(7, fname)

        read(7,*)
        read(7,*) nsubfaults
        read(7,*)

        allocate(subfaults(nsubfaults))

        ! Read in subfaults
        final_rupture_rise_time = 0.d0
        do i=1,nsubfaults
          read(7,*) input_line
          theta = input_line(2)/180.0*4.d0*datan(1.d0)
          subfaults(i)%width = input_line(3)
          subfaults(i)%depth = input_line(4)
          subfaults(i)%slip = input_line(5)
          subfaults(i)%longitude = input_line(9)
          subfaults(i)%rupture_time = input_line(11)
          subfaults(i)%rise_time = input_line(12)

          final_rupture_rise_time = max(final_rupture_rise_time, subfaults(i)%rupture_time + subfaults(i)%rise_time)
        end do

        xp1 = subfaults(1)%longitude*LAT2METER
        yp1 = -subfaults(1)%depth

        xp2 = subfaults(nsubfaults)%longitude*LAT2METER + dcos(theta)*subfaults(nsubfaults)%width
        yp2 = -subfaults(nsubfaults)%depth - dsin(theta)*subfaults(nsubfaults)%width

        center(1) = 0.5d0*(xp1 + xp2)
        center(2) = 0.5d0*(yp1 + yp2)
        total_width = dsqrt((xp2-xp1)**2 + (yp2-yp1)**2)

        xcb(1) = center(1) - 0.5d0*total_width
        xcb(2) = center(1) + 0.5d0*total_width

        ! Determine subfault location in computational domain
        print *, "*******WARNING*******"
        print *, "Subfaults are assumed to be specified in top-center coords"
        print *, "*********************"
        do i=1,nsubfaults
          subfaults(i)%xcb(1) = xcb(1) + dsqrt( &
                (subfaults(i)%longitude*LAT2METER - xp1)**2 + &
                (-subfaults(i)%depth - yp1)**2 )
          subfaults(i)%xcb(2) = subfaults(i)%xcb(1) + subfaults(i)%width
        end do

        close(7)

    end subroutine load_fault

end module fault_module

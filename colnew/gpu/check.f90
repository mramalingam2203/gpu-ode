module testmod
    implicit none
    type A_t
        integer :: u = 0, v = 0
    endtype
    type B_t
        integer :: x = 0, y = 0
        type(A_t) :: a
    endtype

    real a

contains
        a = 1

    subroutine setval ( elem, val )
        integer :: elem(:), val    !! assumed-shape array
        elem(:) = val
    endsubroutine

    subroutine setval_2D ( elem, val )
        integer :: elem(:,:), val    !! assumed-shape array
        elem(:,:) = val
    endsubroutine

    subroutine setval_explicit ( elem, n, val )
        integer :: n, elem( n ), val  !! explicit-shape array
        elem( 1:n ) = val
    endsubroutine
end module

program main
    use testmod
    implicit none
    type(B_t), target :: b( 2 ), bmat( 2, 2 )

    !! Pass 1D component arrays.
    call setval ( b(:)% x,     1 )
    call setval ( b(:)% y,     2 )
    call setval ( b(:)% a% u,  3 )
    call setval ( b(:)% a% v,  4 )

    print *, "b( : )% x    = ", b( : )% x
    print *, "b( : )% y    = ", b( : )% y
    print *, "b( : )% a% u = ", b( : )% a% u
    print *, "b( : )% a% v = ", b( : )% a% v
    print *, "b(1) = ", b(1)
    print *, "b(2) = ", b(2)

    !! Pass a 2D component array.
    call setval_2D ( bmat(:,:)% x, 50 )

    print *, "bmat(:,:)% x = ", bmat(:,:)% x

    !! Pass 1D component array sections.
    call setval ( bmat(:, 1)% x, 1 )
    call setval ( bmat(:, 2)% x, 2 )
    call setval ( bmat(2, :)% x, 5 )

    print *, "bmat(:,:)% x = ", bmat(:,:)% x

    !! Pass a 2D component array to an explicit-shape dummy array
    !! (in this case, copy-in/copy-out may occur).
    call setval_explicit ( bmat(:,:)% x, size(bmat), 100 )

    print *, "bmat(:,:)% x = ", bmat(:,:)% x
endprogram
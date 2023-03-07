# Euclidean metric

## Environment

System:   Ubuntu 22.04
Complier: GNU Fortran 11.3.0
Profiler: GNU prof 2.38

## Method.1

```fortran
dx = x(:, i) - x(:, j)
r  = norm2(dx)
```

## Method.2

```fortran
dx = x(:, i) - x(:, j)
dr = sum(dx*dx)
r  = sqrt(dr)
```

## Method.3

```fortran
dx(1) = x(1, i) - x(1, j)
dr = dx(1)*dx(1)
do d = 2, dim
    dx(d) = x(d, i) - x(d, j)
    dr = dr + dx(d)*dx(d)
end do
r = sqrt(dr)
```

## Code snippet from link list search subroutine

```fortran
...
do
    if (j > i) then
        !!! Method.1
        dx = x(:, i) - x(:, j)
        r  = norm2(dx)
        !!! Method.2
        dx = x(:, i) - x(:, j)
        dr = sum(dx*dx)
        r = sqrt(dr)

        !!! Method.3
        dx(1) = x(1, i) - x(1, j)
        dr = dx(1)*dx(1)
        do d = 2, dim
            dx(d) = x(d, i) - x(d, j)
            dr = dr + dx(d)*dx(d)
        end do
        r = sqrt(dr)
        if (r < hsml(i) * scale_k) then
            if (niac < max_interaction) then
                !!! Neighboring neighborList list, and total interaction number
                !!! and the interaction number for each particle
                niac = niac + 1
                neighborList(niac, 1) = i
                neighborList(niac, 2) = j
                neighborNum(i) = neighborNum(i) + 1
                neighborNum(j) = neighborNum(j) + 1
                !!! Kernel and derivations of kernel
                call kernel(r, dx, hsml(i), w(niac), dwdx(:, niac))
            else
                call print_error(niac, &
                    'Number of interaction pairs exceeds maxium('// &
                    to_string(max_interaction)//')', &
                    "runtime")
                error stop "at subroutine link_list_search"
            end if
        end if
        j = celldata(j)
    else
        exit
    end if
end do
...
```

## Profiling for Method.1

```text
...
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 93.37    125.70   125.70    10000     0.01     0.01  __nnps_m_MOD_link_list_search
...
-----------------------------------------------
              125.70    1.60   10000/10000       __nnps_m_MOD_search_particles [4]
[5]     94.6  125.70    1.60   10000         __nnps_m_MOD_link_list_search [5]
                0.31    1.08 113413519/151813519     __kernel_m_MOD_kernel [7]
                0.19    0.00 19200000/19200000     __nnps_m_MOD_grid_geometry [12]
                0.01    0.00      40/10040       __in_force_m_MOD_in_force [6]
-----------------------------------------------
...
```

## Profiling for Method.2

```test
...
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 90.31     86.32    86.32    10000     0.01     0.01  __nnps_m_MOD_link_list_search
...
-----------------------------------------------
               86.32    1.63   10000/10000       __nnps_m_MOD_search_particles [4]
[5]     92.0   86.32    1.63   10000         __nnps_m_MOD_link_list_search [5]
                0.41    1.06 113413519/151813519     __kernel_m_MOD_kernel [7]
                0.14    0.00 19200000/19200000     __nnps_m_MOD_grid_geometry [13]
                0.01    0.00      40/10040       __in_force_m_MOD_in_force [6]
-----------------------------------------------
...
```

## Profiling for Method.3

```text
...
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 85.50     54.96    54.96    10000     0.01     0.01  __nnps_m_MOD_link_list_search
...
-----------------------------------------------
               54.96    1.71   10000/10000       __nnps_m_MOD_search_particles [4]
[5]     88.2   54.96    1.71   10000         __nnps_m_MOD_link_list_search [5]
                0.35    1.11 113413519/151813519     __kernel_m_MOD_kernel [7]
                0.23    0.00 19200000/19200000     __nnps_m_MOD_grid_geometry [11]
                0.01    0.00      40/10040       __in_force_m_MOD_in_force [6]
-----------------------------------------------
...
```

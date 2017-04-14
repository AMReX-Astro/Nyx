













module amrex_omp_module

  implicit none

  integer, external :: omp_get_num_threads
  integer, external :: omp_get_max_threads
  integer, external :: omp_get_thread_num
  logical, external :: omp_in_parallel

end module amrex_omp_module



program toml_test
   implicit none

   call test_access_file()
   call test_sph_file()

end program toml_test

subroutine test_access_file()
   use tomlf
   use parse_toml_m, only: in_path, out_path
   implicit none
   type(toml_table), allocatable :: access_file
   type(toml_table), pointer :: subtable
   integer :: flag

   ! call system('pwd')

   open(newunit=flag, file='../../access.toml', status='old')
   call toml_parse(access_file, flag)
   close(flag)

   call get_value(access_file, 'access', subtable)
   call get_value(subtable, 'Path', in_path, './data')
   out_path = in_path//"/"//"output"
   write(*,*) in_path
   write(*,*) out_path

   nullify(subtable)

end subroutine test_access_file


subroutine test_sph_file()
   use tomlf
   use parse_toml_m, only: in_path, project_name
   implicit none
   type(toml_table), allocatable :: sph_file
   type(toml_table), pointer :: subtable
   integer :: flag


   open(newunit=flag, file="../."//in_path//'/sph.toml', status='old')
   call toml_parse(sph_file, flag)
   close(flag)

   call get_value(sph_file, 'name', project_name, 'untitled')  !! 默认项目名
   call get_value(sph_file, 'parameter', subtable)

   ! call get_value(subtable, 'kpair', kpair, 20)


   write(*,*) project_name
   ! write(*,*) kpair

   nullify(subtable)

end subroutine test_sph_file

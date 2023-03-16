program toml_test
   implicit none

   call test_access_file()
   call test_sph_file()

end program toml_test

subroutine test_access_file()
   use tomlf
   use ctrl_dict, only: Project
   implicit none
   type(toml_table), allocatable :: access_file
   type(toml_table), pointer :: subtable
   integer :: flag

   ! call system('pwd')

   open(newunit=flag, file='../../access.toml', status='old')
   call toml_parse(access_file, flag)
   close(flag)

   call get_value(access_file, 'access', subtable)
   call get_value(subtable, 'Path', Project%in_path, './data')
   Project%out_path = Project%in_path//"/"//"output"
   write(*,*) Project%in_path
   write(*,*) Project%out_path

   nullify(subtable)

end subroutine test_access_file


subroutine test_sph_file()
   use tomlf
   use ctrl_dict, only: Project
   implicit none
   type(toml_table), allocatable :: sph_file
   type(toml_table), pointer :: subtable
   integer :: flag


   open(newunit=flag, file="../."//Project%in_path//'/sph.toml', status='old')
   call toml_parse(sph_file, flag)
   close(flag)

   call get_value(sph_file, 'name', Project%project_name, 'untitled')  !! 默认项目名
   call get_value(sph_file, 'parameter', subtable)

   ! call get_value(subtable, 'kpair', kpair, 20)


   write(*,*) Project%project_name
   ! write(*,*) kpair

   nullify(subtable)

end subroutine test_sph_file

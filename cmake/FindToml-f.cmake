find_path( Toml-f_INCLUDE_DIR
    NAMES tomlf.mod
    PATHS ~/Repository/include
)

find_library( Toml-f_LIBRARY
    NAMES toml-f
    PATHS ~/Repository/lib
)


set(LOCAL_TOMLF ON CACHE BOOL "User specifies whether to use the local package")
if( Toml-f_INCLUDE_DIR AND Toml-f_LIBRARY AND LOCAL_TOMLF)
    set( Toml-f_FOUND TRUE )
endif()
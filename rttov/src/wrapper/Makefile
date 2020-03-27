
TOPDIR=../..
SUBDIR=wrapper

FFLAGS_ARCH_rttov_wrapper_f2py           = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_wrapper_handle         = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_wrapper_transfer       = $(FFLAGS_ARCH)

F2PYFLAGS_ARCH_rttov_wrapper_f2py       = $(F2PYFLAGS_ARCH)

F90_INTF=rttov_wrapper_f2py.F90

INTF=$(INCDIR)/rttov_wrapper_f2py.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/rttov_wrapper_f2py$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/rttov_wrapper_handle$(EXT_OBJ) $(OBJDIR)/rttov_wrapper_transfer$(EXT_OBJ)

OBJBIN=

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_wrapper$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_brdf_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_other$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_mw_scatt$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_emis_atlas$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_parallel$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(LIBDIR)/rttov_wrapper_f2py$(EXT_SO)

MODS=$(MODDIR)/rttov_wrapper_handle$(EXT_MOD) $(MODDIR)/rttov_wrapper_transfer$(EXT_MOD)

include $(TOPDIR)/build/Makefile.inc

lib: dir intf $(LIB)

bin: lib $(BINS)

dir: 
	@mkdir -p $(TMPDIR)/$(SUBDIR) $(BINDIR) $(OBJDIR) $(INCDIR) $(MODDIR) $(LIBDIR)

intf: 
	@$(MKINTF) $(F90_INTF) $(INTF)
	@$(CPINCH) $(INCDIR)

$(LIB): $(OBJSUB) $(OBJMOD)
	\rm -f $(LIB) $(LIBNOVER)
	$(AR) $(LIB) $(OBJSUB) $(OBJMOD)
	\ln -s $(LIBNAME) $(LIBNOVER)

clean:
	\rm -f $(INTF) $(LIB) $(LIBNOVER) $(BINS) $(OBJSUB) $(OBJMOD) $(OBJBIN) $(MODS)


$(LIBDIR)/rttov_wrapper_f2py$(EXT_SO): rttov_wrapper_f2py.F90 $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_setup_brdf_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_get_brdf.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_emis_ret_terms.interface $(INCDIR)/rttov_scatt_setupindex.interface $(INCDIR)/rttov_scatt.interface $(INCDIR)/rttov_parallel_scatt.interface $(INCDIR)/rttov_alloc_k.interface $(INCDIR)/rttov_init_rad.interface $(INCDIR)/rttov_init_transmission.interface $(INCDIR)/rttov_init_prof.interface $(INCDIR)/rttov_k.interface $(INCDIR)/rttov_parallel_k.interface $(INCDIR)/rttov_init_scatt_prof.interface $(INCDIR)/rttov_scatt_ad.interface $(INCDIR)/rttov_parallel_scatt_ad.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_opts_scatt.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_legcoef_calc.interface $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_wrapper_transfer$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_wrapper_handle$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(LIBS)
	cd $(TMPDIR) && cp -f $(SRCDIR)/f2py_f2cmap .f2py_f2cmap
	cd $(TMPDIR) && $(F2PY) -c $(SRCDIR)/rttov_wrapper_f2py.F90 -m rttov_wrapper_f2py \
	 --f77flags=$(F2PYFLAGS_ARCH_rttov_wrapper_f2py) --f90flags=$(F2PYFLAGS_ARCH_rttov_wrapper_f2py) $(FFLAGS_EXTERN) \
	$(FFLAG_MOD)$(MODDIR) -I$(INCDIR) \
	-L$(LIBDIR) -lrttov$(RTTOV_VERSION)_wrapper -lrttov$(RTTOV_VERSION)_brdf_atlas -lrttov$(RTTOV_VERSION)_other -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_emis_atlas -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_parallel -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(F2PYLDFLAGS_ARCH)
	@cd $(TMPDIR) && $(MVDMOD) rttov_wrapper_f2py.cpython-37m-x86_64-linux-gnu.so $(LIBDIR)/rttov_wrapper_f2py$(EXT_SO)

$(OBJDIR)/rttov_wrapper_f2py$(EXT_OBJ): rttov_wrapper_f2py.F90 $(INCDIR)/rttov_user_options_checkinput.interface $(INCDIR)/rttov_errorreport.interface $(INCDIR)/rttov_setup_emis_atlas.interface $(INCDIR)/rttov_setup_brdf_atlas.interface $(INCDIR)/rttov_get_emis.interface $(INCDIR)/rttov_get_brdf.interface $(INCDIR)/rttov_alloc_direct.interface $(INCDIR)/rttov_direct.interface $(INCDIR)/rttov_parallel_direct.interface $(INCDIR)/rttov_init_opt_param.interface $(INCDIR)/rttov_alloc_scatt_prof.interface $(INCDIR)/rttov_alloc_emis_ret_terms.interface $(INCDIR)/rttov_scatt_setupindex.interface $(INCDIR)/rttov_scatt.interface $(INCDIR)/rttov_parallel_scatt.interface $(INCDIR)/rttov_alloc_k.interface $(INCDIR)/rttov_init_rad.interface $(INCDIR)/rttov_init_transmission.interface $(INCDIR)/rttov_init_prof.interface $(INCDIR)/rttov_k.interface $(INCDIR)/rttov_parallel_k.interface $(INCDIR)/rttov_init_scatt_prof.interface $(INCDIR)/rttov_scatt_ad.interface $(INCDIR)/rttov_parallel_scatt_ad.interface $(INCDIR)/rttov_print_opts.interface $(INCDIR)/rttov_print_opts_scatt.interface $(INCDIR)/rttov_bpr_init.interface $(INCDIR)/rttov_bpr_calc.interface $(INCDIR)/rttov_bpr_dealloc.interface $(INCDIR)/rttov_legcoef_calc.interface $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD) $(MODDIR)/rttov_wrapper_transfer$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_bpr_mod$(EXT_MOD) $(MODDIR)/rttov_htfrtc_interface_mod$(EXT_MOD) $(MODDIR)/rttov_wrapper_handle$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/rttov_scattering_mod$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_wrapper_f2py) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_wrapper_f2py.F90 -o $(OBJDIR)/rttov_wrapper_f2py$(EXT_OBJ)


$(OBJDIR)/rttov_wrapper_handle$(EXT_OBJ): rttov_wrapper_handle.F90 $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_dealloc_scattcoeffs.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_read_scattcoeffs.interface $(INCDIR)/rttov_nullify_scattcoeffs.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_alloc_emis_ret_terms.interface $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_wrapper_handle) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_wrapper_handle.F90 -o $(OBJDIR)/rttov_wrapper_handle$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_wrapper_handle$(EXT_MOD) $(MODDIR)/rttov_wrapper_handle$(EXT_MOD)

$(MODDIR)/rttov_wrapper_handle$(EXT_MOD): rttov_wrapper_handle.F90 $(INCDIR)/rttov_dealloc_coefs.interface $(INCDIR)/rttov_dealloc_scattcoeffs.interface $(INCDIR)/rttov_deallocate_emis_atlas.interface $(INCDIR)/rttov_deallocate_brdf_atlas.interface $(INCDIR)/rttov_read_coefs.interface $(INCDIR)/rttov_read_scattcoeffs.interface $(INCDIR)/rttov_nullify_scattcoeffs.interface $(INCDIR)/rttov_alloc_rad.interface $(INCDIR)/rttov_alloc_transmission.interface $(INCDIR)/rttov_alloc_emis_ret_terms.interface $(MODDIR)/mod_uwiremis_atlas$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/rttov_coef_io_mod$(EXT_MOD) $(MODDIR)/mod_camel_atlas$(EXT_MOD) $(MODDIR)/mod_rttov_emis_atlas$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/mod_rttov_brdf_atlas$(EXT_MOD) $(MODDIR)/mod_brdf_atlas$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/mod_cnrm_mw_atlas$(EXT_MOD) $(MODDIR)/rttov_hdf_mod$(EXT_MOD) $(MODDIR)/mod_mwatlas_m2$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_wrapper_handle) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_wrapper_handle.F90 -o $(OBJDIR)/rttov_wrapper_handle$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_wrapper_handle$(EXT_MOD) $(MODDIR)/rttov_wrapper_handle$(EXT_MOD)


$(OBJDIR)/rttov_wrapper_transfer$(EXT_OBJ): rttov_wrapper_transfer.F90  $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_wrapper_handle$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_wrapper_transfer) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_wrapper_transfer.F90 -o $(OBJDIR)/rttov_wrapper_transfer$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_wrapper_transfer$(EXT_MOD) $(MODDIR)/rttov_wrapper_transfer$(EXT_MOD)

$(MODDIR)/rttov_wrapper_transfer$(EXT_MOD): rttov_wrapper_transfer.F90  $(MODDIR)/rttov_const$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/rttov_wrapper_handle$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_wrapper_transfer) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_wrapper_transfer.F90 -o $(OBJDIR)/rttov_wrapper_transfer$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) rttov_wrapper_transfer$(EXT_MOD) $(MODDIR)/rttov_wrapper_transfer$(EXT_MOD)


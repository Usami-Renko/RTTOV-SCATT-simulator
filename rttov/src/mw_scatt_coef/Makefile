
TOPDIR=../..
SUBDIR=mw_scatt_coef

FFLAGS_ARCH_arts_scat                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_convert_mietable                = $(FFLAGS_ARCH)
FFLAGS_ARCH_density_all                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_gamma_dsd                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_get_dia                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_ice_density                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_liu_dda                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_liu_density                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_load_arts_ssp                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_melting_layer                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_mg_ellips                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_mie_coated_sphere               = $(FFLAGS_ARCH)
FFLAGS_ARCH_mie_one_temp                    = $(FFLAGS_ARCH)
FFLAGS_ARCH_mie_one_wc                      = $(FFLAGS_ARCH)
FFLAGS_ARCH_mie_sphere                      = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_arts                        = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_gamma_dsd                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_mod_mie                         = $(FFLAGS_ARCH)
FFLAGS_ARCH_n0_t                            = $(FFLAGS_ARCH)
FFLAGS_ARCH_perm_ice                        = $(FFLAGS_ARCH)
FFLAGS_ARCH_perm_melt                       = $(FFLAGS_ARCH)
FFLAGS_ARCH_perm_water                      = $(FFLAGS_ARCH)
FFLAGS_ARCH_perm_water_TKC_16               = $(FFLAGS_ARCH)
FFLAGS_ARCH_perm_water_liebe_89             = $(FFLAGS_ARCH)
FFLAGS_ARCH_perm_water_rosenkranz_15        = $(FFLAGS_ARCH)
FFLAGS_ARCH_permittivity                    = $(FFLAGS_ARCH)
FFLAGS_ARCH_predict_mom07                   = $(FFLAGS_ARCH)
FFLAGS_ARCH_predict_psd                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_predict_psd_F07                 = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_ascii2bin_scattcoef       = $(FFLAGS_ARCH)
FFLAGS_ARCH_rttov_scatt_make_coef           = $(FFLAGS_ARCH)
FFLAGS_ARCH_scattering                      = $(FFLAGS_ARCH)
FFLAGS_ARCH_set_spectra                     = $(FFLAGS_ARCH)
FFLAGS_ARCH_vol_fracs                       = $(FFLAGS_ARCH)

CFLAGS_ARCH_scatdb     = $(CFLAGS_ARCH)

F90_INTF=arts_scat.F90 convert_mietable.F90 density_all.F90 gamma_dsd.F90 get_dia.F90 liu_dda.F90 liu_density.F90 load_arts_ssp.F90 melting_layer.F90 mg_ellips.F90 mie_coated_sphere.F90 mie_one_temp.F90 mie_one_wc.F90 mie_sphere.F90 mod_gamma_dsd.F90 n0_t.F90 perm_ice.F90 perm_melt.F90 perm_water.F90 perm_water_TKC_16.F90 perm_water_liebe_89.F90 perm_water_rosenkranz_15.F90 permittivity.F90 predict_mom07.F90 predict_psd.F90 predict_psd_F07.F90 scattering.F90 set_spectra.F90 vol_fracs.F90

INTF=$(INCDIR)/arts_scat.interface $(INCDIR)/convert_mietable.interface $(INCDIR)/density_all.interface $(INCDIR)/gamma_dsd.interface $(INCDIR)/get_dia.interface $(INCDIR)/liu_dda.interface $(INCDIR)/liu_density.interface $(INCDIR)/load_arts_ssp.interface $(INCDIR)/melting_layer.interface $(INCDIR)/mg_ellips.interface $(INCDIR)/mie_coated_sphere.interface $(INCDIR)/mie_one_temp.interface $(INCDIR)/mie_one_wc.interface $(INCDIR)/mie_sphere.interface $(INCDIR)/mod_gamma_dsd.interface $(INCDIR)/n0_t.interface $(INCDIR)/perm_ice.interface $(INCDIR)/perm_melt.interface $(INCDIR)/perm_water.interface $(INCDIR)/perm_water_TKC_16.interface $(INCDIR)/perm_water_liebe_89.interface $(INCDIR)/perm_water_rosenkranz_15.interface $(INCDIR)/permittivity.interface $(INCDIR)/predict_mom07.interface $(INCDIR)/predict_psd.interface $(INCDIR)/predict_psd_F07.interface $(INCDIR)/scattering.interface $(INCDIR)/set_spectra.interface $(INCDIR)/vol_fracs.interface

LIBNAME=librttov$(RTTOV_VERSION)_$(SUBDIR)$(EXT_LIB)
LIB=$(LIBDIR)/$(LIBNAME)
LIBNOVER=$(LIBDIR)/librttov_$(SUBDIR)$(EXT_LIB)

OBJSUB=$(OBJDIR)/arts_scat$(EXT_OBJ) $(OBJDIR)/convert_mietable$(EXT_OBJ) $(OBJDIR)/density_all$(EXT_OBJ) $(OBJDIR)/gamma_dsd$(EXT_OBJ) $(OBJDIR)/get_dia$(EXT_OBJ) $(OBJDIR)/ice_density$(EXT_OBJ) $(OBJDIR)/liu_dda$(EXT_OBJ) $(OBJDIR)/liu_density$(EXT_OBJ) $(OBJDIR)/load_arts_ssp$(EXT_OBJ) $(OBJDIR)/melting_layer$(EXT_OBJ) $(OBJDIR)/mg_ellips$(EXT_OBJ) $(OBJDIR)/mie_coated_sphere$(EXT_OBJ) $(OBJDIR)/mie_one_temp$(EXT_OBJ) $(OBJDIR)/mie_one_wc$(EXT_OBJ) $(OBJDIR)/mie_sphere$(EXT_OBJ) $(OBJDIR)/mod_gamma_dsd$(EXT_OBJ) $(OBJDIR)/n0_t$(EXT_OBJ) $(OBJDIR)/perm_ice$(EXT_OBJ) $(OBJDIR)/perm_melt$(EXT_OBJ) $(OBJDIR)/perm_water$(EXT_OBJ) $(OBJDIR)/perm_water_TKC_16$(EXT_OBJ) $(OBJDIR)/perm_water_liebe_89$(EXT_OBJ) $(OBJDIR)/perm_water_rosenkranz_15$(EXT_OBJ) $(OBJDIR)/permittivity$(EXT_OBJ) $(OBJDIR)/predict_mom07$(EXT_OBJ) $(OBJDIR)/predict_psd$(EXT_OBJ) $(OBJDIR)/predict_psd_F07$(EXT_OBJ) $(OBJDIR)/scattering$(EXT_OBJ) $(OBJDIR)/set_spectra$(EXT_OBJ) $(OBJDIR)/vol_fracs$(EXT_OBJ) $(OBJDIR)/scatdb$(EXT_OBJ) 

OBJMOD=$(OBJDIR)/mod_arts$(EXT_OBJ) $(OBJDIR)/mod_mie$(EXT_OBJ)

OBJBIN=$(OBJDIR)/rttov_ascii2bin_scattcoef$(EXT_OBJ) $(OBJDIR)/rttov_scatt_make_coef$(EXT_OBJ)

LIBS=$(LIBDIR)/librttov$(RTTOV_VERSION)_mw_scatt_coef$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_mw_scatt$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

BINS=$(BINDIR)/rttov_ascii2bin_scattcoef$(EXT_BIN) $(BINDIR)/rttov_scatt_make_coef$(EXT_BIN)

MODS=$(MODDIR)/mod_arts$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)

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


$(BINDIR)/rttov_ascii2bin_scattcoef$(EXT_BIN): $(OBJDIR)/rttov_ascii2bin_scattcoef$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_ascii2bin_scattcoef$(EXT_BIN) $(OBJDIR)/rttov_ascii2bin_scattcoef$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_mw_scatt_coef -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(BINDIR)/rttov_scatt_make_coef$(EXT_BIN): $(OBJDIR)/rttov_scatt_make_coef$(EXT_OBJ) $(LIBS)
	cd $(TMPDIR) && $(FC) -o $(BINDIR)/rttov_scatt_make_coef$(EXT_BIN) $(OBJDIR)/rttov_scatt_make_coef$(EXT_OBJ) -L$(LIBDIR) -lrttov$(RTTOV_VERSION)_mw_scatt_coef -lrttov$(RTTOV_VERSION)_mw_scatt -lrttov$(RTTOV_VERSION)_coef_io -lrttov$(RTTOV_VERSION)_hdf -lrttov$(RTTOV_VERSION)_main $(LDFLAGS_EXTERN) $(LDFLAGS_ARCH)

$(OBJDIR)/scatdb$(EXT_OBJ): scatdb.c
	cd $(TMPDIR) && $(CC) -I$(INCDIR) $(CFLAGS_EXTERN) $(CFLAGS_ARCH_scatdb) -c $(SRCDIR)/scatdb.c -o $(OBJDIR)/scatdb$(EXT_OBJ)


$(OBJDIR)/arts_scat$(EXT_OBJ): arts_scat.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_arts_scat) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/arts_scat.F90 -o $(OBJDIR)/arts_scat$(EXT_OBJ)


$(OBJDIR)/convert_mietable$(EXT_OBJ): convert_mietable.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_convert_mietable) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/convert_mietable.F90 -o $(OBJDIR)/convert_mietable$(EXT_OBJ)


$(OBJDIR)/density_all$(EXT_OBJ): density_all.F90 $(INCDIR)/liu_density.interface $(MODDIR)/mod_arts$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_density_all) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/density_all.F90 -o $(OBJDIR)/density_all$(EXT_OBJ)


$(OBJDIR)/gamma_dsd$(EXT_OBJ): gamma_dsd.F90 $(INCDIR)/n0_t.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_gamma_dsd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/gamma_dsd.F90 -o $(OBJDIR)/gamma_dsd$(EXT_OBJ)


$(OBJDIR)/get_dia$(EXT_OBJ): get_dia.F90  $(MODDIR)/mod_arts$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_get_dia) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/get_dia.F90 -o $(OBJDIR)/get_dia$(EXT_OBJ)


$(OBJDIR)/ice_density$(EXT_OBJ): ice_density.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_ice_density) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/ice_density.F90 -o $(OBJDIR)/ice_density$(EXT_OBJ)


$(OBJDIR)/liu_dda$(EXT_OBJ): liu_dda.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_liu_dda) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/liu_dda.F90 -o $(OBJDIR)/liu_dda$(EXT_OBJ)


$(OBJDIR)/liu_density$(EXT_OBJ): liu_density.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_liu_density) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/liu_density.F90 -o $(OBJDIR)/liu_density$(EXT_OBJ)


$(OBJDIR)/load_arts_ssp$(EXT_OBJ): load_arts_ssp.F90  $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_load_arts_ssp) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/load_arts_ssp.F90 -o $(OBJDIR)/load_arts_ssp$(EXT_OBJ)


$(OBJDIR)/melting_layer$(EXT_OBJ): melting_layer.F90 $(INCDIR)/perm_melt.interface $(INCDIR)/mie_coated_sphere.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_melting_layer) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/melting_layer.F90 -o $(OBJDIR)/melting_layer$(EXT_OBJ)


$(OBJDIR)/mg_ellips$(EXT_OBJ): mg_ellips.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mg_ellips) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mg_ellips.F90 -o $(OBJDIR)/mg_ellips$(EXT_OBJ)


$(OBJDIR)/mie_coated_sphere$(EXT_OBJ): mie_coated_sphere.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mie_coated_sphere) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mie_coated_sphere.F90 -o $(OBJDIR)/mie_coated_sphere$(EXT_OBJ)


$(OBJDIR)/mie_one_temp$(EXT_OBJ): mie_one_temp.F90 $(INCDIR)/get_dia.interface $(INCDIR)/scattering.interface $(INCDIR)/set_spectra.interface $(INCDIR)/mie_one_wc.interface $(MODDIR)/mod_arts$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mie_one_temp) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mie_one_temp.F90 -o $(OBJDIR)/mie_one_temp$(EXT_OBJ)


$(OBJDIR)/mie_one_wc$(EXT_OBJ): mie_one_wc.F90 $(INCDIR)/permittivity.interface $(INCDIR)/melting_layer.interface $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mie_one_wc) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mie_one_wc.F90 -o $(OBJDIR)/mie_one_wc$(EXT_OBJ)


$(OBJDIR)/mie_sphere$(EXT_OBJ): mie_sphere.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mie_sphere) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mie_sphere.F90 -o $(OBJDIR)/mie_sphere$(EXT_OBJ)


$(OBJDIR)/mod_gamma_dsd$(EXT_OBJ): mod_gamma_dsd.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_gamma_dsd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_gamma_dsd.F90 -o $(OBJDIR)/mod_gamma_dsd$(EXT_OBJ)


$(OBJDIR)/n0_t$(EXT_OBJ): n0_t.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_n0_t) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/n0_t.F90 -o $(OBJDIR)/n0_t$(EXT_OBJ)


$(OBJDIR)/perm_ice$(EXT_OBJ): perm_ice.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_perm_ice) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/perm_ice.F90 -o $(OBJDIR)/perm_ice$(EXT_OBJ)


$(OBJDIR)/perm_melt$(EXT_OBJ): perm_melt.F90 $(INCDIR)/vol_fracs.interface $(INCDIR)/mg_ellips.interface $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_perm_melt) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/perm_melt.F90 -o $(OBJDIR)/perm_melt$(EXT_OBJ)


$(OBJDIR)/perm_water$(EXT_OBJ): perm_water.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_perm_water) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/perm_water.F90 -o $(OBJDIR)/perm_water$(EXT_OBJ)


$(OBJDIR)/perm_water_TKC_16$(EXT_OBJ): perm_water_TKC_16.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_perm_water_TKC_16) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/perm_water_TKC_16.F90 -o $(OBJDIR)/perm_water_TKC_16$(EXT_OBJ)


$(OBJDIR)/perm_water_liebe_89$(EXT_OBJ): perm_water_liebe_89.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_perm_water_liebe_89) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/perm_water_liebe_89.F90 -o $(OBJDIR)/perm_water_liebe_89$(EXT_OBJ)


$(OBJDIR)/perm_water_rosenkranz_15$(EXT_OBJ): perm_water_rosenkranz_15.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_perm_water_rosenkranz_15) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/perm_water_rosenkranz_15.F90 -o $(OBJDIR)/perm_water_rosenkranz_15$(EXT_OBJ)


$(OBJDIR)/permittivity$(EXT_OBJ): permittivity.F90 $(INCDIR)/perm_water_liebe_89.interface $(INCDIR)/perm_water_rosenkranz_15.interface $(INCDIR)/perm_water_TKC_16.interface $(INCDIR)/perm_ice.interface $(INCDIR)/mg_ellips.interface $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_permittivity) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/permittivity.F90 -o $(OBJDIR)/permittivity$(EXT_OBJ)


$(OBJDIR)/predict_mom07$(EXT_OBJ): predict_mom07.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_predict_mom07) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/predict_mom07.F90 -o $(OBJDIR)/predict_mom07$(EXT_OBJ)


$(OBJDIR)/predict_psd$(EXT_OBJ): predict_psd.F90 $(INCDIR)/density_all.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_predict_psd) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/predict_psd.F90 -o $(OBJDIR)/predict_psd$(EXT_OBJ)


$(OBJDIR)/predict_psd_F07$(EXT_OBJ): predict_psd_F07.F90 $(INCDIR)/predict_mom07.interface $(INCDIR)/density_all.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_predict_psd_F07) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/predict_psd_F07.F90 -o $(OBJDIR)/predict_psd_F07$(EXT_OBJ)


$(OBJDIR)/scattering$(EXT_OBJ): scattering.F90 $(INCDIR)/mie_sphere.interface $(INCDIR)/permittivity.interface $(INCDIR)/liu_dda.interface $(INCDIR)/density_all.interface $(INCDIR)/arts_scat.interface $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_scattering) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/scattering.F90 -o $(OBJDIR)/scattering$(EXT_OBJ)


$(OBJDIR)/set_spectra$(EXT_OBJ): set_spectra.F90 $(INCDIR)/predict_psd.interface $(INCDIR)/predict_psd_F07.interface $(INCDIR)/gamma_dsd.interface $(INCDIR)/mod_gamma_dsd.interface $(INCDIR)/density_all.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_set_spectra) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/set_spectra.F90 -o $(OBJDIR)/set_spectra$(EXT_OBJ)


$(OBJDIR)/vol_fracs$(EXT_OBJ): vol_fracs.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_vol_fracs) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/vol_fracs.F90 -o $(OBJDIR)/vol_fracs$(EXT_OBJ)


$(OBJDIR)/rttov_ascii2bin_scattcoef$(EXT_OBJ): rttov_ascii2bin_scattcoef.F90 $(INCDIR)/rttov_opencoeff.interface $(INCDIR)/rttov_read_scattcoeffs.interface $(INCDIR)/rttov_errorreport.interface $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/rttov_global$(EXT_MOD) $(MODDIR)/rttov_unix_env$(EXT_MOD) $(MODDIR)/rttov_getoptions$(EXT_MOD) $(MODDIR)/rttov_types$(EXT_MOD) $(MODDIR)/yomhook$(EXT_MOD) $(MODDIR)/rttov_const$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_ascii2bin_scattcoef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_ascii2bin_scattcoef.F90 -o $(OBJDIR)/rttov_ascii2bin_scattcoef$(EXT_OBJ)


$(OBJDIR)/rttov_scatt_make_coef$(EXT_OBJ): rttov_scatt_make_coef.F90 $(INCDIR)/convert_mietable.interface $(INCDIR)/mie_one_temp.interface $(INCDIR)/permittivity.interface $(INCDIR)/liu_dda.interface $(INCDIR)/load_arts_ssp.interface $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_rttov_scatt_make_coef) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/rttov_scatt_make_coef.F90 -o $(OBJDIR)/rttov_scatt_make_coef$(EXT_OBJ)


$(OBJDIR)/mod_arts$(EXT_OBJ): mod_arts.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_arts) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_arts.F90 -o $(OBJDIR)/mod_arts$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_arts$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)

$(MODDIR)/mod_arts$(EXT_MOD): mod_arts.F90  $(MODDIR)/mod_mie$(EXT_MOD) $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_arts) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_arts.F90 -o $(OBJDIR)/mod_arts$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_arts$(EXT_MOD) $(MODDIR)/mod_arts$(EXT_MOD)


$(OBJDIR)/mod_mie$(EXT_OBJ): mod_mie.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_mie) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_mie.F90 -o $(OBJDIR)/mod_mie$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_mie$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)

$(MODDIR)/mod_mie$(EXT_MOD): mod_mie.F90  $(MODDIR)/parkind1$(EXT_MOD)
	cd $(TMPDIR) && $(FC) -I$(INCDIR) $(FFLAGS_EXTERN) $(FFLAGS_ARCH_mod_mie) $(FFLAG_MOD)$(MODDIR) -c $(SRCDIR)/mod_mie.F90 -o $(OBJDIR)/mod_mie$(EXT_OBJ)
	@cd $(TMPDIR) && $(MVDMOD) mod_mie$(EXT_MOD) $(MODDIR)/mod_mie$(EXT_MOD)


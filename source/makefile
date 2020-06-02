FC = ifort
FCOPTS = -O3 
LD = ifort
LDOPTS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#LDOPTS = -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXENAME = HF_run
OBJS =HFB.o read_input.o dimm.o deallocate_all.o initialize_basis.o occupations.o kinetic.o interaction.o transition.o hfb_iteration.o transf_interaction.o HFB_energy.o hp_field.o hn_field.o Dp_field.o Dn_field.o hY_field.o diagonalization.o make_sp_levels.o make_densities.o check_densities.o R_val.o gauss_int.o order.o TDA.o NYTDA.o pnTDAhh.o pnTDAhp.o pnTDAph.o pnTDApp.o ppTDApp.o nnTDApp.o TDAtz.o TDAtzpp.o interaction_NY.o make_qsp_levels.o MBPT_energy.o transit_calc.o transit_calc_spur.o phonon_density_calc.o diag_subspace.o spur_vec.o spur_vectz.o geom.o technical.o math.o geometric.o 

$(EXENAME): $(OBJS)
	$(LD) $(LDOPTS) -o ../$(EXENAME) $(OBJS)

HFB.o: HFB.f technical.o math.o
	$(FC) $(FCOPTS) -c HFB.f

read_input.o: read_input.f
	$(FC) $(FCOPTS) -c read_input.f

dimm.o: dimm.f
	$(FC) $(FCOPTS) -c dimm.f

deallocate_all.o: deallocate_all.f technical.o
	$(FC) $(FCOPTS) -c deallocate_all.f

interaction.o: interaction.f technical.o geom.o
	$(FC) $(FCOPTS) -c interaction.f

transition.o: transition.f technical.o geom.o
	$(FC) $(FCOPTS) -c transition.f

hfb_interaction.o: hfb_interaction.f technical.o geom.o
	$(FC) $(FCOPTS) -c hfb_interaction.f

initialize_basis.o: initialize_basis.f technical.o
	$(FC) $(FCOPTS) -c initialize_basis.f
	
occupations.o: occupations.f technical.o
	$(FC) $(FCOPTS) -c occupations.f

kinetic.o: kinetic.f technical.o
	$(FC) $(FCOPTS) -c kinetic.f

transf_interaction.o: transf_interaction.f technical.o
	$(FC) $(FCOPTS) -c transf_interaction.f

HFB_energy.o: HFB_energy.f technical.o
	$(FC) $(FCOPTS) -c HFB_energy.f

hp_field.o: hp_field.f technical.o
	$(FC) $(FCOPTS) -c hp_field.f

hn_field.o: hn_field.f technical.o
	$(FC) $(FCOPTS) -c hn_field.f

Dp_field.o: Dp_field.f technical.o
	$(FC) $(FCOPTS) -c Dp_field.f

Dn_field.o: Dn_field.f technical.o
	$(FC) $(FCOPTS) -c Dn_field.f

hY_field.o: hY_field.f technical.o
	$(FC) $(FCOPTS) -c hY_field.f
	
make_sp_levels.o: make_sp_levels.f technical.o
	$(FC) $(FCOPTS) -c make_sp_levels.f

make_densities.o: make_densities.f technical.o
	$(FC) $(FCOPTS) -c make_densities.f

check_densities.o: check_densities.f technical.o
	$(FC) $(FCOPTS) -c check_densities.f

R_val.o: R_val.f technical.o
	$(FC) $(FCOPTS) -c R_val.f

gauss_int.o: gauss_int.f technical.o
	$(FC) $(FCOPTS) -c gauss_int.f

order.o: order.f technical.o
	$(FC) $(FCOPTS) -c order.f

diagonalization.o: diagonalization.f
	$(FC) $(FCOPTS) -c diagonalization.f

TDA.o: TDA.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c TDA.f

NYTDA.o: NYTDA.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c NYTDA.f

pnTDApp.o: pnTDApp.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c pnTDApp.f

pnTDAhh.o: pnTDAhh.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c pnTDAhh.f

pnTDAph.o: pnTDAph.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c pnTDAph.f

pnTDAhp.o: pnTDAhp.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c pnTDAhp.f

ppTDApp.o: ppTDApp.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c ppTDApp.f

nnTDApp.o: nnTDApp.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c nnTDApp.f

TDAtz.o: TDAtz.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c TDAtz.f

TDAtzpp.o: TDAtzpp.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c TDAtzpp.f

interaction_NY.o: interaction_NY.f technical.o geometric.o
	$(FC) $(FCOPTS) -c interaction_NY.f

make_qsp_levels.o: make_qsp_levels.f technical.o
	$(FC) $(FCOPTS) -c make_qsp_levels.f

MBPT_energy.o: MBPT_energy.f technical.o geom.o
	$(FC) $(FCOPTS) -c MBPT_energy.f

transit_calc.o: transit_calc.f technical.o
	$(FC) $(FCOPTS) -c transit_calc.f

transit_calc_spur.o: transit_calc_spur.f technical.o
	$(FC) $(FCOPTS) -c transit_calc_spur.f

phonon_density_calc.o: phonon_density_calc.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c phonon_density_calc.f

diag_subspace.o: diag_subspace.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c diag_subspace.f

spur_vec.o: spur_vec.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c spur_vec.f

spur_vectz.o: spur_vectz.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c spur_vectz.f

technical.o: technical.f
	$(FC) $(FCOPTS) -c technical.f

math.o: math.f
	$(FC) $(FCOPTS) -c math.f

geom.o: geom.f
	$(FC) $(FCOPTS) -c geom.f

geometric.o: geometric.f90
	$(FC) $(FCOPTS) -c geometric.f90

#.f90.f.o:
#	$(FC) $(FCOPTS) -c $<

clean:
	rm -f *.o *.mod ../$(EXENAME)

#clobber:
#	rm -f $(EXENAME)

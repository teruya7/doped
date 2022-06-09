from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.defects.generators import VacancyGenerator, SubstitutionGenerator, VoronoiInterstitialGenerator
from pymatgen.analysis.defects.core import Interstitial
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from doped.pycdt.core.defectsmaker import DefectChargerSemiconductor, \
     SimpleChargeGenerator, DefectChargerInsulator, \
     DefectChargerUserCustom, DefectChargerIonic


class ChargedDefectsSurfaces(object):
    """
    A class to generate charged defective surface slabs for use in first
    principles supercell formalism. The standard defects such as vacancies, 
    antisites, substitutions and interstitals are generated. The method will not work for 
    slabs without inversion symmetry.
    """
    def __init__(self, structure,  max_min_oxi=None, substitutions=None,
                 oxi_states=None, vacancies_flag=True, antisites_flag=True, 
                 interstitials_flag=False, inter_elements=None, 
                 struct_type='semiconductor', dev=5):
        """
        Args:
            structure (Slab):
                the supercell expansion of the slab, must have oxidation 
                states added. 
            max_min_oxi (dict):
                The minimal and maximum oxidation state of each element as a
                dict. For instance {"O":(-2,0)}. If not given, the oxi-states
                of pymatgen are considered.
            substitutions (dict):
                The allowed substitutions of elements as a dict. If not given,
                intrinsic defects are computed. If given, intrinsic (e.g.,
                anti-sites) and extrinsic are considered explicitly specified.
                Example: {"Co":["Zn","Mn"]} means Co sites can be substituted
                by Mn or Zn.
            oxi_states (dict):
                The oxidation state of the elements in the compound e.g.
                {"Fe":2,"O":-2}. Required for checking the slabs are non-polar
            vacancies_flag (bool):
                If False, don't generate vacancies.
            antisites_flag (bool):
                If False, don't generate antisites.
            interstitals_flag (bool):
                If True, generate interstitials
            inter_elements (list): A list of strings containing symbols that are to be considered for interstitial sites. The default (None)
                triggers self-interstitial generation,
                given that include_interstitials is True.
            struct_type (string):
                Options are 'semiconductor' and 'insulator'. If semiconductor
                is selected, charge states based on database of semiconductors
                is used to assign defect charges. For insulators, defect
                charges are conservatively assigned.
            dev (float): The maximum distance in angstroms from the surface of the slab to the defect, default is 5 Å
        """
        max_min_oxi = max_min_oxi if max_min_oxi is not None else {}
        substitutions = substitutions if substitutions is not None else {}
        oxi_states = oxi_states if oxi_states is not None else {}
        inter_elements = inter_elements if inter_elements is not None else []

        # check all inputs are legit
        if interstitials_flag and inter_elements:
            for elem_str in inter_elements:
                if not Element.is_valid_symbol(elem_str):
                    raise ValueError("invalid interstitial element"
                            " \"{}\"".format(elem_str))

        if not isinstance(structure, Slab): 
            raise TypeError('Structure must be pymatgen.core.surface.Slab') 

        self.defects = []
        self.substitutions = {}
        self.struct_type = struct_type
        for key, val in substitutions.items():
            self.substitutions[key] = val

        self.struct = structure

        # find where the coords of zmin and zmax 
        self.zmax = max(self.struct.cart_coords[:,2])
        self.zmin = min(self.struct.cart_coords[:,2])

        if self.struct_type == 'semiconductor':
            self.defect_charger = DefectChargerSemiconductor(self.struct,
                                                             min_max_oxi=max_min_oxi, 
                                                             oxi_states=oxi_states)
        elif self.struct_type == 'insulator':
            self.defect_charger = DefectChargerInsulator(self.struct)
        elif self.struct_type == 'manual':
            self.defect_charger = DefectChargerUserCustom(self.struct,
                                                          oxi_states=oxi_states)
        elif self.struct_type == 'ionic':
            self.defect_charger = DefectChargerIonic(self.struct)
        else:
            raise NotImplementedError

        self.defects = {}
        sc = self.struct.copy()
        sc_scale = (1,1,1)
        self.defects['bulk'] = {'name': 'bulk',
                'supercell': {'size': sc_scale, 'structure': sc}}

        vac_defs = []
        as_defs = []
        sub_defs = []
        ints_defs = []

        if vacancies_flag: 
            print('generating vacancies')
            vg = VacancyGenerator(self.struct)
            for i, vac in enumerate(vg): 
                symbol = vac.site.specie.symbol
                charges_vac = self.defect_charger.get_charges('vacancy', symbol)
                if vac.site.coords[2] - self.zmin < dev or self.zmax - vac.site.coords[2] < dev:
                    vac_sc = self.struct.copy() 
                    # index in the host for one of the two sites, remove on both sides
                    index = self.struct.sites.index(vac.site)
                    vac_sc.add_oxidation_state_by_element(oxi_states)
                    vac_sc.symmetrically_remove_atoms([index])
                    vac_sc.add_oxidation_state_by_element(oxi_states)

                    struct_for_defect_site = Structure(vac.bulk_structure.copy().lattice,
                                                    [vac.site.specie],
                                                    [vac.site.frac_coords],
                                                    to_unit_cell=True, coords_are_cartesian=False)
                    vac_sc_site = struct_for_defect_site[0]

                    # check for invesion & polarity 
                    if vac_sc.is_symmetric() and not vac_sc.is_polar(): 
                        for c in SimpleChargeGenerator(vac):
                            vacancy = { 
                                'name': f'vac_{i+1}_{symbol}',
                                'site_specie': symbol,
                                'site_multiplicity': vac.multiplicity,
                                'bulk_supercell_site': vac_sc_site,
                                'unique_site': vac.site, 
                                'charges': charges_vac,
                                'supercell': {'size': sc_scale, 'structure': vac_sc},
                                'index': index, 
                                'defect_type': 'vacancy', 
                                'Possible_KV_charge': c.charge
                            }
                            vac_defs.append(vacancy)

        if antisites_flag: 
            print('generating antisites')
            for as_specie in list(self.struct.types_of_specie):
                SG = SubstitutionGenerator(self.struct, as_specie)
                for i, sub in enumerate(SG):
                    if sub.site.coords[2] - self.zmin < dev or self.zmax - sub.site.coords[2] < dev:
                        # create a trivial defect structure to find where 
                        # supercell transformation moves the defect
                        struct_for_defect_site = Structure(sub.bulk_structure.copy().lattice,
                                                        [sub.site.specie],
                                                        [sub.site.frac_coords],
                                                        to_unit_cell=True, 
                                                        coords_are_cartesian=False)
                        as_sc_site = struct_for_defect_site[0]

                        poss_deflist = sorted(sub.bulk_structure.get_sites_in_sphere(sub.site.coords, 
                                            0.01, include_index=True), key=lambda x: x[1])

                        defindex = poss_deflist[0][2]
                        as_site = sub.bulk_structure[defindex] # antisite site in bulk structure
                        vac_symbol = as_site.specie.symbol # symbol of the removed atom 
                        as_symbol = as_specie.symbol # symbol of the added atom
                        charges_as = self.defect_charger.get_charges('antisite', vac_symbol, as_symbol)
                        
                        # index in the host for one of the two sites, remove on both sides
                        as_sc = self.struct.copy() 
                        index = self.struct.sites.index(as_site)
                        as_sc.add_oxidation_state_by_element(oxi_states)
                        as_sc.symmetrically_remove_atoms([index])
                        
                        as_sc.symmetrically_add_atom(as_symbol, as_site.frac_coords)
                        as_sc.add_oxidation_state_by_element(oxi_states)
                        
            
                        # check if symmetric - if symmetric should also be non-polar but can't check 
                        # for polarity because the fractional coordinate the method gives 
                        # is negative, so the method pymatgen uses would give a polar 
                        # slab because of how it calculates the dipole 
                        # (same applies for interstitial generation)
                        if as_sc.is_symmetric():
                            for c in SimpleChargeGenerator(sub):
                                as_defs.append({
                                            'name': "as_{}_{}_on_{}".format(
                                                i+1, as_symbol, vac_symbol),
                                            'unique_site': as_site,
                                            'bulk_supercell_site': as_sc_site,
                                            'defect_type': 'antisite',
                                            'site_specie': vac_symbol,
                                            'substituting_specie': as_symbol,
                                            'site_multiplicity': sub.multiplicity,
                                            'charges': charges_as,
                                            'Possible_KV_Charge': c.charge,
                                            'supercell': {'size': sc_scale,'structure': as_sc},
                                            })
        
        for vac_symbol, subspecie_list in self.substitutions.items():
            for subspecie_symbol in subspecie_list:
                SG = SubstitutionGenerator(self.struct, subspecie_symbol)
                for i, sub in enumerate(SG):
                    if sub.site.coords[2] - self.zmin < dev or self.zmax - sub.site.coords[2] < dev:
                        # get bulk_site (non sc)
                        poss_deflist = sorted(sub.bulk_structure.get_sites_in_sphere(sub.site.coords, 
                        0.1, include_index=True), key=lambda x: x[1])
                        if not len(poss_deflist):
                            raise ValueError("Could not find substitution site inside" 
                            " bulk structure for {}?".format(sub.name))
                        defindex = poss_deflist[0][2]
                        sub_site = self.struct[defindex] # sub site in bulk structure
                        this_vac_symbol = sub_site.specie.symbol # symbol of removed site
                        sub_symbol = sub.site.specie.symbol # symbol of added site

                        if (sub_symbol != subspecie_symbol) or (this_vac_symbol != vac_symbol):
                            continue
                        else:
                            sub_sc = self.struct.copy() 
                            index = self.struct.sites.index(as_site)
                            sub_sc.symmetrically_remove_atoms([index])
                            sub_sc.symmetrically_add_atom(as_symbol, as_site.frac_coords)
                            sub_sc.add_oxidation_state_by_element(oxi_states)

                            # create a trivial defect structure to find where supercell
                            # transformation moves the defect
                            struct_for_defect_site = Structure(sub.bulk_structure.copy().lattice,
                                                                [sub.site.specie],
                                                                [sub.site.frac_coords],
                                                                to_unit_cell=True, 
                                                                coords_are_cartesian=False)
                            struct_for_defect_site.make_supercell(sc_scale)
                            sub_sc_site = struct_for_defect_site[0]
                            
                            # get possible charges 
                            charges_sub = self.defect_charger.get_charges(
                                    'substitution', vac_symbol, subspecie_symbol)
                            
                            # check for symmetry 
                            if sub_sc.is_symmetric():         
                                for c in SimpleChargeGenerator(sub):
                                    sub_defs.append({
                                        'name': "sub_{}_{}_on_{}".format(
                                            i+1, subspecie_symbol, vac_symbol),
                                        'unique_site': sub_site,
                                        'bulk_supercell_site': sub_sc_site,
                                        'defect_type':'substitution',
                                        'site_specie':vac_symbol,
                                        'substitution_specie':subspecie_symbol,
                                        'site_multiplicity': sub.multiplicity,
                                        'supercell':{'size': sc_scale,'structure': sub_sc},
                                        'charges': charges_sub,
                                        'Possible_KV_Charge': c.charge})

        if interstitials_flag: 
            if not inter_elements:
                inter_elements = [elem.symbol for elem in \
                        self.struct.composition.elements]
            
            print('searching for Voronoi interstitial sites, this can take a while')
            # get the sites for the interstitals for an imaginary interstitial only once
            IG = VoronoiInterstitialGenerator(self.struct, 'O') 
            sites = []
            for ints in IG: 
                # check that the interstitital site is within dev range of the surface
                if ints.site.coords[2] - self.zmin < dev or self.zmax - ints.site.coords[2] < dev:
                    sites.append(ints)

            # iterate over the inter elements to get the correct periodic sites
            inter_sites = []
            for el in inter_elements: 
                for s in sites: 
                    d = PeriodicSite(el, s.site.frac_coords, s.bulk_structure.copy().lattice)
                    inter_sites.append(d)
            
            # make the defect structure  
            for i, intersite in enumerate(inter_sites):
                el = intersite.specie
                name = "inter_{}_{}".format(i+1, el)

                if intersite.lattice != self.struct.lattice:
                    err_msg = "Lattice matching error occurs between provided interstitial and the bulk structure."
                    raise ValueError(err_msg)
                else:
                    ints = Interstitial(self.struct, intersite) 

                # create a trivial defect structure to find where supercell transformation moves the defect site
                struct_for_defect_site = Structure(ints.bulk_structure.copy().lattice,[ints.site.specie], [ints.site.frac_coords],to_unit_cell=True, coords_are_cartesian=False)

                struct_for_defect_site.make_supercell(sc_scale)
                ints_site = struct_for_defect_site[0]

                ints_sc = self.struct.copy() 

                # index in the host for one of the two sites, remove on both sides
                ints_sc.symmetrically_add_atom(el, ints_site.frac_coords)
                ints_sc.add_oxidation_state_by_element(oxi_states)
                charges_ints = self.defect_charger.get_charges(
                        'interstitial', el)

                if ints_sc.is_symmetric():
                    for c in SimpleChargeGenerator(ints):
                        ints_defs.append({
                                    'name': name,
                                    'unique_site': ints.site,
                                    'bulk_supercell_site': ints_site,
                                    'defect_type': 'interstitial',
                                    'site_specie': ints.site.specie.symbol,
                                    'site_multiplicity': ints.multiplicity,
                                    'charges': charges_ints,
                                    'Possible_KV_Charge': c.charge,
                                    'supercell': {'size': sc_scale,'structure': ints_sc},
                                    })


        self.defects['vacancies'] = vac_defs
        self.defects['substitutions'] = sub_defs
        self.defects['substitutions'] += as_defs
        self.defects['interstitials'] = ints_defs

        print("\nNumber of jobs created:")
        tottmp=0
        for j in self.defects.keys():
            if j=='bulk':
                print("    bulk = 1")
                tottmp += 1
            else:
                print("    {}:".format(j))
                for lis in self.defects[j]:
                    print("        {} = {} with site multiplicity {}".format(lis['name'],
                                                   len(lis['charges']), lis['site_multiplicity']))
                    tottmp += len(lis['charges'])
        print("Total (non dielectric) jobs created = {}\n".format(tottmp))
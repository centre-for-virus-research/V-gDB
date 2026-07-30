from django.urls import path
from api.routes import sequences
from api.routes import alignments
from api.routes import mutations
from api.routes import search
from api.routes import phylogeny
from api.routes import lineage
from api.routes import taxonomy
from api.routes import analysis
from api.routes import polymorphisms

urlpatterns = [
    
    # SEQUENCES
    path('sequences/', sequences.get_sequences, name='get_sequences'), #DONE
    path('sequence/<str:primary_accession>', sequences.get_sequence, name='get_sequence'), #done
    path('sequence/reference/<str:primary_accession>', sequences.get_reference_sequence, name='get_reference_sequence'), #done
    path('sequences/global/', sequences.get_global_distribution_of_sequences, name='get_global_distribution_of_sequences'),
    path('sequences/download_sequences_meta_data/', sequences.download_sequences_meta_data, name='download_sequences_meta_data'),
    path('sequences/download_sequences/', sequences.download_sequences, name='download_sequences'),

    # PHYLOGENY
    path("phylogeny/trees/", phylogeny.get_trees, name='get_trees'),

    # POLYMORPHISMS
    path('polymorphisms/', polymorphisms.get_polymorphisms, name='get_polymorphisms'), 
    path('polymorphism/<str:id>', polymorphisms.get_polymorphism, name='get_polymorphism'), 

    # TAXANOMY
    path('taxonomy/phylum', taxonomy.get_phylum, name='get_phylum'), 
    path('taxonomy/class', taxonomy.get_class, name='get_class'), 
    path('taxonomy/order_category', taxonomy.get_order, name='get_order'), 
    path('taxonomy/family', taxonomy.get_family, name='get_family'),
    path('taxonomy/genus', taxonomy.get_genus, name='get_genus'),
    path('taxonomy/species', taxonomy.get_species, name='get_species'),



    # MUTATIONS 
    path('mutations/host_adaptation', mutations.get_host_adaptations, name='get_host_adaptations'),


    path('lineages/', lineage.get_lineage, name='get_lineage'),



    path('alignments/download', alignments.download_alignments, name='download_alignment'),




    path('analysis/clade_assignment/', analysis.run_phylogenetic_clade_assignment_analysis, name='run_phylogenetic_clade_assignment_analysis'),
    path('analysis/drug_analysis/', analysis.run_drug_analysis, name='run_drug_analysis'),
    
    
    # FILTERS
    path('filters/search_isolate_ids/<str:query>', search.search_isolate_ids, name='search_isolate_ids'),
    path('filters/search_pubmed_ids/<str:query>', search.search_pubmed_ids, name='search_pubmed_ids'),
    path('filters/search_hosts/', search.search_hosts, name='search_hosts'),
    path('filters/search_primary_accession_ids/<str:query>', search.search_primary_accession_ids, name='search_primary_accession_ids'),
    path('filters/search_country/', search.search_country, name='search_country'),

    path('filters/search_segment/', search.search_segments, name='search_segments'),


    path('filters/search_m49_intermediate/', search.search_m49_intermediate, name='search_m49_intermediate'),
    path('filters/search_m49_region/', search.search_m49_region, name='search_m49_region'),
    path('filters/search_m49_sub_region/', search.search_m49_sub_region, name='search_m49_sub_region'),



    path('filters/search_protein_name/', search.search_protein_name, name='search_protein_name'),
    path('filters/search_drug/', search.search_drug, name='search_drug'),

    
]
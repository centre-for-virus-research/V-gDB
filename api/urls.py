from django.urls import path
from api.routes import sequences
from api.routes import alignments
from api.routes import versions
from api.routes import mutations
from api.routes import tasks
from api.routes import search
from api.routes import phylogeny
from api.routes import strains
from api.routes import lineage
from api.routes import taxonomy
from api.routes import analysis

urlpatterns = [
    
    # SEQUENCES
    # API paths to get all of the data with various filters
    path('sequences/', sequences.get_sequences, name='get_sequences'), #DONE
    path('sequence/<str:primary_accession>', sequences.get_sequence, name='get_sequence'), #done
    path('sequence/reference/<str:primary_accession>', sequences.get_reference_sequence, name='get_reference_sequence'), #done
    


    path('strains/', strains.get_strains, name='get_strains'), 
    path('strain/<path:strain_id>', strains.get_strain, name='get_strain'), #done


    path('lineages/', lineage.get_lineage, name='get_lineage'),



    path('alignments/download', alignments.download_alignments, name='download_alignment'),
    path('sequences/download_sequences_meta_data/', sequences.download_sequences_meta_data, name='download_sequences_meta_data'),
    path('sequences/download_sequences/', sequences.download_sequences, name='download_sequences'),



    # Taxonomy
    path('taxonomy/phylum', taxonomy.get_phylum, name='get_phylum'), 
    path('taxonomy/class', taxonomy.get_class, name='get_class'), 
    path('taxonomy/order_category', taxonomy.get_order, name='get_order'), 
    path('taxonomy/family', taxonomy.get_family, name='get_family'),
    path('taxonomy/genus', taxonomy.get_genus, name='get_genus'),
    path('taxonomy/species', taxonomy.get_species, name='get_species'),
    


    # NOT FINISHED STUFF
    path('sequences/global/', sequences.get_global_distribution_of_sequences, name='get_global_distribution_of_sequences'),
    path("phylogeny/tree/", phylogeny.get_tree, name='get_tree'),

    path('adaptive_mutations/', mutations.get_adaptive_mutations, name='get_adaptive_mutations'),

    path('analysis/mutations/', mutations.get_mutations, name='get_mutations'),



    path('analysis/clade_assignment/', analysis.run_phylogenetic_clade_assignment_analysis, name='run_phylogenetic_clade_assignment_analysis'),

    
    # API paths to get data for primary accession
    
    
    # MUTATIONS 
    path('adaptive_mutations/', mutations.get_adaptive_mutations_chart, name='get_adaptive_mutations_chart'),

    path('adaptive_mutations_chart/<str:segment>', mutations.get_adaptive_mutations_chart, name='get_adaptive_mutations_chart'),

    

    
    
   


    path('get_host_species/', sequences.get_host_species, name='get_host_species'),

    # ALIGNMENTS - DONE (documentation API)
    path('alignments/download_alignments/', alignments.download_alignments, name='download_alignments'),


    # MUTATIONS - DONE (documentation API)
    
    path('mutations/get_mutation_regions_and_codons', mutations.get_mutation_regions_and_codons, name='get_mutations_regions_and_codons'),
    


    
    # TASKS
    path('tasks/run_sequence_alignment/', tasks.run_sequence_alignment, name='run_sequence_alignment'),
    path('tasks/get_blast_results/<str:job_id>', tasks.get_blast_results, name='get_blast_results'),
    
    
    path('tasks/get_alignment_results/<str:job_id>', tasks.get_alignment_results, name='get_alignment_results'),
    
    
    
    path('tasks/get_job_logs/<str:job_id>', tasks.get_job_logs, name='get_job_logs'),







    # STATISTICS - DONE (documentation API)
    

    



    # FILTERS
    path('filters/search_isolate_ids/<str:query>', search.search_isolate_ids, name='search_isolate_ids'),
    path('filters/search_pubmed_ids/<str:query>', search.search_pubmed_ids, name='search_pubmed_ids'),
    path('filters/search_hosts/', search.search_hosts, name='search_hosts'),
    path('filters/search_primary_accession_ids/<str:query>', search.search_primary_accession_ids, name='search_primary_accession_ids'),
    path('filters/search_country/', search.search_country, name='search_country'),


    path('filters/search_m49_intermediate/', search.search_m49_intermediate, name='search_m49_intermediate'),
    path('filters/search_m49_region/', search.search_m49_region, name='search_m49_region'),
    path('filters/search_m49_sub_region/', search.search_m49_sub_region, name='search_m49_sub_region'),



    path('filters/search_region/', search.search_region, name='search_region'),
    path('filters/search_phylum/', search.search_phylum, name='search_phylum'),



    path('check_db_connection', versions.check_db_connection, name='check_db_connection'),

    
    
    path('advanced_filter/<str:query>', sequences.advanced_filter, name='advanced_filter'),



    path('get_vgt_version/', versions.get_vgt_version, name='get_vgt_version'),
    path('get_meta_data_columns/', versions.get_meta_data_columns, name='get_meta_data_columns')
    
]
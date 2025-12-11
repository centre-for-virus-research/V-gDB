from django.urls import path
from api.routes import sequences
from api.routes import alignments
from api.routes import versions
from api.routes import mutations
from api.routes import statistics
from api.routes import features
from api.routes import genes
from api.routes import tasks
from api.routes import search
from api.routes import phylogeny
# MAIN API endpoints! 
# Virus Web Resource GUI will 

urlpatterns = [
    
    # SEQUENCES
    # API paths to get all of the data with various filters
    path('sequences/', sequences.get_sequences, name='get_sequences'), 
    path('sequences/metadata/', sequences.get_sequences_meta_data, name='get_sequences_meta_data'),
    path('sequences/alignment/', sequences.get_sequences_alignment, name='get_sequences_alignment'),
    path('sequences/reference/metadata/', sequences.get_reference_sequences_meta_data, name='get_reference_sequences_meta_data'),


    path('sequences/metadata/map/', sequences.get_map_metadata, name='get_map_metadata'),
    
    # API paths to get data for primary accession
    path('sequence/metadata/<str:primary_accession>', sequences.get_sequence_meta_data, name='get_sequence_meta_data'), #done
    path('sequence/alignment/<str:primary_accession>', sequences.get_sequence_alignment, name='get_sequence_alignment'),  #done
    path('sequence/reference/<str:primary_accession>', sequences.get_reference_sequence, name='get_reference_sequence'), #done




    path('sequences/strains/', sequences.get_strains, name='get_strains'), 
    path('sequence/strain/<path:isolate>', sequences.get_strain, name='get_strain'), #done


    path("phylogeny/tree/", phylogeny.get_tree, name='get_tree'),








    # MUTATIONS 
    path('adaptive_mutations/', mutations.get_adaptive_mutations, name='get_adaptive_mutations'),

    path('analysis/mutations/', mutations.get_mutations, name='get_mutations'),

    path('alignments/get_reference_sequence/<str:primary_accession>', sequences.get_reference_sequence, name='get_reference_sequence'),
    
    
    path('sequences/download_sequences_meta_data/', sequences.download_sequences_meta_data, name='download_sequences_meta_data'),


    path('get_host_species/', sequences.get_host_species, name='get_host_species'),

    # ALIGNMENTS - DONE (documentation API)
    path('alignments/download_alignments/', alignments.download_alignments, name='download_alignments'),


    # MUTATIONS - DONE (documentation API)
    path('mutations/get_mutations/', mutations.get_mutations, name='get_mutations'),
    path('mutations/get_mutation_regions_and_codons', mutations.get_mutation_regions_and_codons, name='get_mutations_regions_and_codons'),
    


    
    # TASKS
    path('tasks/run_sequence_alignment/', tasks.run_sequence_alignment, name='run_sequence_alignment'),
    path('tasks/get_blast_results/<str:job_id>', tasks.get_blast_results, name='get_blast_results'),
    
    
    path('tasks/get_alignment_results/<str:job_id>', tasks.get_alignment_results, name='get_alignment_results'),
    
    
    
    path('tasks/get_job_logs/<str:job_id>', tasks.get_job_logs, name='get_job_logs'),







    # STATISTICS - DONE (documentation API)
    path('statistics/get_global_distribution_of_sequences/', statistics.get_global_distribution_of_sequences, name='get_global_distribution_of_sequences'),
    path('statistics/get_statistics/', statistics.get_statistics, name='get_statistics'),

    





    # Features
    path('features/get_features/', features.get_features, name='get_features'),

    # Genes
    path('genes/get_genes_tree/', genes.get_genes_tree, name='get_genes_tree'),

    path('features/get_clades_tree/', features.get_clades_tree, name='get_clades_tree'),


    # FILTERS
    path('filters/search_isolate_ids/<str:query>', search.search_isolate_ids, name='search_isolate_ids'),
    path('filters/search_pubmed_ids/<str:query>', search.search_pubmed_ids, name='search_pubmed_ids'),
    path('filters/search_hosts/<str:query>', search.search_hosts, name='search_hosts'),
    path('filters/search_primary_accession_ids/<str:query>', search.search_primary_accession_ids, name='search_primary_accession_ids'),
    path('filters/search_country/<str:query>', search.search_country, name='search_country'),
    path('filters/search_region/', search.search_region, name='search_region'),



    path('check_db_connection', versions.check_db_connection, name='check_db_connection'),

    
    
    path('advanced_filter/<str:query>', sequences.advanced_filter, name='advanced_filter'),



    path('get_vgt_version/', versions.get_vgt_version, name='get_vgt_version'),
    path('get_meta_data_columns/', versions.get_meta_data_columns, name='get_meta_data_columns')
    
]
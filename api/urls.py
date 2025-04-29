from django.urls import path
from api.routes import sequences
from api.routes import alignments
from api.routes import versions
from api.routes import mutations
from api.routes import statistics
from api.routes import features
from api.routes import genes
from api.routes import tasks
from api.routes import filters

# MAIN API endpoints! 
# Virus Web Resource GUI will 

urlpatterns = [


    # SEQUENCES
    path('sequences/get_sequences_meta_data/', sequences.get_sequences_meta_data, name='get_sequences_meta_data'),
    path('sequences/get_sequence_meta_data/<str:primary_accession>', sequences.get_sequence_meta_data, name='get_sequence_meta_data'),
    path('alignments/get_reference_sequences_meta_data/', sequences.get_reference_sequences_meta_data, name='get_reference_sequences_meta_data'),
    path('alignments/get_reference_sequence/<str:primary_accession>', sequences.get_reference_sequence, name='get_reference_sequence'),
    
    
    path('sequences/get_sequences_meta_data_by_filters/', sequences.get_sequences_meta_data_by_filters, name='get_sequences_meta_data_by_filters'),
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


    # FILTERS
    path('filters/search_isolate_ids/<str:query>', filters.search_isolate_ids, name='search_isolate_ids'),
    path('filters/search_pubmed_ids/<str:query>', filters.search_pubmed_ids, name='search_pubmed_ids'),
    path('filters/search_hosts/<str:query>', filters.search_hosts, name='search_hosts'),
    path('filters/search_primary_accession_ids/<str:query>', filters.search_primary_accession_ids, name='search_primary_accession_ids'),



    path('check_db_connection', versions.check_db_connection, name='check_db_connection'),

    
    
    path('advanced_filter/<str:query>', sequences.advanced_filter, name='advanced_filter'),



    path('get_vgt_version/', versions.get_vgt_version, name='get_vgt_version'),
    path('get_meta_data_columns/', versions.get_meta_data_columns, name='get_meta_data_columns')
    
]
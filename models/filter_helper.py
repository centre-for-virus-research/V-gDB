# from models import sequences_helpers as sh
from django.db import connections

TAXONOMY_FILTERS = ['phylum', 'class', 'order_category', 'family', 'genus', 'species', 'host']
REGION_FILTERS = ['country', 'm49_region_id', 'm49_sub_region_id', 'm49_intermediate_region_id', 'm49_code']
GENOME_FILTERS = ['full_genome', 'nucleoprotein', 'phosphoprotein', 'm2_protein', 'glycoprotein', 'l_protein', 'coreprotein', 'envelope_protein_E1', 'envelope_protein_E2', 'protein_p7', 'NS2', 'NS3', 'NS4A', 'NS4B', 'NS5B']
CLADE_FILTERS = ['EPA_major_clade', 'EPA_minor_clade']
STANDARD_FILTERS = ['primary_accession', 'isolate', 'exclusion_status', 'accession_type', 'country_validated']
COMPARISON_FILTERS = {
                'length_lower': ('real_length', '<='),
                'length_upper': ('real_length', '>='),
                'collection_year_lower': ('collection_year', '<='),
                'collection_year_upper': ('collection_year', '>='),
                'creation_year_lower': ('create_date', '>='),
                'creation_year_upper': ('create_date', '<=')
            }


FILTER_GROUPS = {
                    **{k: "comparison" for k in COMPARISON_FILTERS},
                    **{k: "taxonomy" for k in TAXONOMY_FILTERS},
                    **{k: "region" for k in REGION_FILTERS},
                    **{k: "genome" for k in GENOME_FILTERS},
                    **{k: "clade" for k in CLADE_FILTERS},
                    **{k: "standard" for k in STANDARD_FILTERS},
                }

PROTEIN_FILTERS_MAP = {
                        'nucleoprotein': 'nucleoprotein N',
                        'phosphoprotein': 'phosphoprotein M1',
                        'm2_protein': 'M2 protein',
                        'glycoprotein': 'transmembrane glycoprotein G',
                        'l_protein': 'L protein',
                        'coreprotein': 'core protein',
                        'envelope_protein_E1': 'envelope protein E1',
                        'envelope_protein_E2': 'envelope protein E2',
                        'protein_p7': 'protein p7',
                        'NS2': 'nonstructural protein NS2',
                        'NS3': 'protease/helicase protein NS3',
                        'NS4A': 'nonstructural protein NS4A',
                        'NS4B': 'nonstructural protein NS4B',
                        'NS5B': 'RNA-dependent RNA polymerase NS5B'
                    }

class FilterHelper:

    def __init__(self, filters, database):

        self.filters = filters
        self.database = database
        self.filter_keys = []
        self.where_clauses, self.params = [], []
        self.taxonomy_clauses, self.taxonomy_params = [], []
        self.region_clauses, self.region_params = [], []
        self.genome_coverage_clauses, self.genome_coverage_params = [], []

        self.handlers = {
                        "comparison": self._handle_comparison,
                        "taxonomy": self._handle_taxonomy,
                        "region": self._handle_region,
                        "genome": self._handle_genome,
                        "clade": self._handle_clade,
                        "standard": self._handle_standard,
                    }

    def add_filters(self):

        self.filter_keys = self.filters.keys()
        for key, value in self.filters.items():
            key, exclude = self._parse_exclude(key)

            group = FILTER_GROUPS.get(key)
            if group:
                self.handlers[group](key, value, exclude)

        self._build_taxonomy()
        self._build_region()
        self._build_genome_coverage()

        where_str = ' AND '.join(self.where_clauses)

        print(where_str, self.params)

        return where_str, self.params
    
    def add_filters_host_mutations(self):
        clauses, params, columns = [], [], []
        for key, value in self.filters.items():
            columns.append(f"hl.{key} as host")
            if isinstance(value, list):
                placeholders = ', '.join(['%s'] * len(value))
                clauses.append(f"(hl.{key} IN ({placeholders}))")
                params.extend(value)
            else:
                clauses.append(f"(hl.{key} = %s)")
                params.append(value)

        where_str = ' AND '.join(clauses)
        columns_str = ', '.join(columns)
        return where_str, params, columns_str


    def _parse_exclude(self, key):
        if key.startswith("exclude_"):
            return key[8:], True
        return key, False
    
    def _get_exclusion_comparison(self, key):
        comparison = 'IN'
        if (key in self.filter_keys):
            comparison = 'NOT IN'
        return comparison

    
    def _handle_taxonomy(self, key, value, exclude):
        clause, param = self._add_standard_filters(key, value, exclude)
        self.taxonomy_clauses.append(clause)
        self.taxonomy_params.extend(param)

    def _handle_comparison(self, key, value, exclude):
        col, op = COMPARISON_FILTERS[key]
        self.where_clauses.append(f"{col} {op} %s")
        self.params.append(int(value))

    def _handle_region(self, key, value, exclude):
        if key=='country':
            key='display_name'
        clause, param = self._add_standard_filters(key, value, exclude)
        self.region_clauses.append(clause)
        self.region_params.extend(param)

    def _handle_genome(self, key, value, exclude):
        if key == 'full_genome':
            clause, param = self._add_standard_filters('calculated_genome_coverage', value, exclude)
            self.where_clauses.append(clause)
            self.params.extend(param)
        else:

            product = PROTEIN_FILTERS_MAP.get(key)
            if product:
                clause, param = self._add_genome_coverage_filters(product, value)
                self.genome_coverage_clauses = clause
                self.genome_coverage_params.extend(param)
            print(self.genome_coverage_clauses)

    def _handle_clade(self, key, value, exclude):
        if self.filters.get("exclude_clades"):
            exclude = True
            clause, param = self._add_exclude_clade_filters(key, self.filters.get("EPA_major_clade"), self.filters.get("EPA_minor_clade"), exclude)
        else:
            clause, param = self._add_standard_filters(key, value, exclude)
        
        self.where_clauses.append(clause)
        self.params.extend(param)

    def _handle_standard(self, key, value, exclude):
        clause, param = self._add_standard_filters(key, value, exclude)
        self.where_clauses.append(clause)
        self.params.extend(param)

    def _build_taxonomy(self):
        if self.taxonomy_clauses:
            comparison = self._get_exclusion_comparison("exclude_taxa")
            if ("host" in self.filter_keys):
                clause, param = self._add_taxonomy_host_filters(self.filters["host"])
                self.params.extend(param)
            else:
                clause = self._add_taxonomy_filters(self.taxonomy_clauses, comparison)
                self.params.extend(self.taxonomy_params)
                
            self.where_clauses.append(clause)

    def _build_region(self):

        if self.region_clauses:

            comparison = self._get_exclusion_comparison("exclude_region")
            
            clause = self._add_region_filters(self.region_clauses, comparison)
            self.params.extend(self.region_params)
            self.where_clauses.append(clause)

    def _build_genome_coverage(self):

        if self.genome_coverage_clauses:
            comparison = self._get_exclusion_comparison("exclude_genome")
            clause = self._add_genome_filter(self.genome_coverage_clauses, comparison)
            self.params.extend(self.genome_coverage_params)
            self.where_clauses.append(clause)

    def _get_comparison(self, value, exclude):
        if isinstance(value, list):
            comparison = "IN"
            if (exclude):
                comparison = "NOT IN"
        else:
            comparison = "="
            if (exclude):
                comparison = "!="
        return comparison

    def _add_standard_filters(self, key, value, exclude):

        params, where_clauses, placeholders = [], [], []
        comparison = self._get_comparison(value, exclude)

        if isinstance(value, list):
            placeholders = ', '.join(['%s'] * len(value))
            where_clauses = f"{key} {comparison} ({placeholders})"
            params.extend(value)
        else:
            where_clauses = f"{key} {comparison} %s"
            params.append(value)

        return where_clauses, params
    
    def _add_genome_coverage_filters(self, product, coverage):

        params, where_clauses = [], []
        where_clauses.append("product = %s")
        where_clauses.append("genome_coverage >= %s")
        params.append(product)
        params.append(coverage)

        return where_clauses, params
    
    def _add_genome_filter(self, clauses, comparison):
        print(clauses)
        where_str = ' AND '.join(clauses)
        print(where_str)
        sql = f"""
                primary_accession {comparison} (
                    SELECT accession
                    FROM features
                    WHERE {where_str}
                )
                """
        return sql
    
    def _add_exclude_clade_filters(self, key, major, minor=None, exclude=True):

        params, where_clauses, placeholders_major, placeholders_minor = [], [], [], []
        # comparison = get_comparison(value, exclude)

        if isinstance(major, list):
            placeholders_major = ', '.join(['%s'] * len(major))
            where_clauses_major = f"EPA_major_clade IN ({placeholders_major})"
            params.extend(major)
        else:
            where_clauses_major = f"EPA_major_clade = %s"
            params.append(major)

        if minor:
            if isinstance(minor, list):
                placeholders_minor = ', '.join(['%s'] * len(minor))
                where_clauses_minor = f"EPA_minor_clade IN ({placeholders_minor})"
                params.extend(minor)
            else:
                where_clauses_minor = f"EPA_minor_clade = %s"
                params.append(minor)

        if minor:
            where_clauses = f"NOT ({where_clauses_major} AND {where_clauses_minor})"
        else:
            where_clauses = f"NOT ({where_clauses_major})"
            

        return where_clauses, params
    
    def _add_taxonomy_host_filters(self, value):

        params = []
        common_clause, common_param = self._add_standard_filters("common_name", value, False)
        common_sql = f""" SELECT taxa_id FROM host_taxa WHERE {common_clause} """
        
        params.extend(common_param)

        recursive_sql = f""" WITH RECURSIVE taxa_tree(id) AS (
                        SELECT lineage_taxa_id
                        FROM host_lineage_lookup
                        WHERE lineage_taxa_id IN ({common_sql})

                    UNION

                        SELECT hll.desc_taxa_id
                        FROM host_lineage_lookup hll
                        JOIN taxa_tree tt
                        ON hll.lineage_taxa_id = tt.id 
                    )
                SELECT DISTINCT id FROM taxa_tree;
                """
        
        with connections[self.database].cursor() as cursor:
            cursor.execute(recursive_sql, params)
            rows = cursor.fetchall()
            ids = [row[0] for row in rows]

        # return recursive_sql, params
        clause, param = self._add_standard_filters("host_taxa_id", ids, False)

        return clause, param

    def _add_taxonomy_filters(self, clauses, comparison):
        taxa_where_str = ' AND '.join(clauses)
        taxa_sql = f""" 
                        host_taxa_id {comparison} (
                                                    SELECT taxa_id
                                                    FROM host_lineage
                                                    WHERE {taxa_where_str}
                                                  )
            """
        return taxa_sql
    
    def _add_region_filters(self, clauses, comparison):

        region_where_str = ' AND '.join(clauses)
        region_sql = f"""
                        country_validated {comparison} (
                                                            SELECT m49_code
                                                            FROM m49_country 
                                                            WHERE {region_where_str}
                                                        )
                    """
        return region_sql
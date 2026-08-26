API - Sequences
===============

sequences
---------

This API request will return a list of all the metadata and sequence
information for the corresponding filters. If no filters are added, it
will return a list of every sequence available in the database.

Usage
~~~~~

.. code-block:: bash

   python script.py --taxid <TAXID> [options]

Arguments
~~~~~~~~~

**isolate**

**(Optional)** Isolate ID of the sequence.

Example: ``11292``


**primary_accession**

**(Optional)** NCBI accession of the sequence.

Example: ``AB009601``


**pubmed_id**

**(Optional)** PubMed ID number that is used as a reference for the sequence.

Example: ``11270607``


**exclusion_status**

**(Optional)** Some sequences are automatically excluded due to various
reasons. To exclude these, set ``exclusion_status`` to ``0``. If you only
want excluded sequences, set ``exclusion_status`` to ``1``.

Example: ``0, 1``


Genome Coverage
^^^^^^^^^^^^^^^

**(Optional)** This filter allows you to find sequences based on their
genome coverage.

The following genome coverage filters are available:

**Greater Than**

Finds all sequences that have a genome coverage greater than the
specified percentage.

Example: ``90``


**Greater Than or Equal To**

Finds all sequences that have a genome coverage greater than or equal to
the specified percentage.

Example: ``90``


**Less Than**

Finds all sequences that have a genome coverage less than the specified
percentage.

Example: ``50``


**Less Than or Equal To**

Finds all sequences that have a genome coverage less than or equal to
the specified percentage.

Example: ``50``


**Equal To**

Finds all sequences that have a genome coverage equal to the specified
percentage.

Example: ``100``


Example
~~~~~~~

.. code-block:: bash

   python script.py --taxid 11292 --email myname@domain.com --batch_size 200 --update
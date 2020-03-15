DROP TABLE IF EXISTS results;
DROP TABLE IF EXISTS tree;

CREATE TABLE results (
  id TEXT PRIMARY KEY,
    sample_name TEXT NOT NULL,
    status TEXT DEFAULT "processing",
    project_id TEXT,
    result TEXT,
    created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    user_id TEXT,
    branch TEXT
);

CREATE TABLE tree (
    created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    newick TEXT
);

CREATE TABLE tree_data (
    accession TEXT,
    collection_data TIMESTAMP,
    country TEXT
)

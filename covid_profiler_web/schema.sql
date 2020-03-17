DROP TABLE IF EXISTS results;
DROP TABLE IF EXISTS tree;
DROP TABLE IF EXISTS tree_data;
DROP TABLE IF EXISTS mutations;

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
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    newick TEXT
);

CREATE TABLE tree_data (
    id TEXT,
    collection_date TEXT,
    country TEXT
)

CREATE TABLE mutations (
    position INT,
    origins INT,
    samples TEXT
)

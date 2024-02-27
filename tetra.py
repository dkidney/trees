from ete3 import Tree
import pandas as pd

df = pd.DataFrame([
    ('Amniota', 'Root', 100),
    ('Synapsida', 'Amniota', 320),
    ('Sauropsida', 'Amniota', 75),
    # ('Eureptilia', 'Sauropsida', 10),
    # ('Parareptilia', 'Sauropsida', 33),
    # ('Diapsids', 'Eureptilia', 10),
    ('non-Archosauria', 'Sauropsida', 245),
    ('Archosauria', 'Sauropsida', 95),
    # ('Pseudosuchia', 'Archosauria', 10),
    # ('Avemetatarsalia', 'Archosauria', 10),
    ('Birds', 'Archosauria', 150),
    ('non-Birds', 'Archosauria', 150),
], columns=['node', 'parent', 'branch.length'])

print(df)

parent_child_list = df[['parent', 'node', 'branch.length']].apply(tuple, axis=1).tolist()
tree = Tree.from_parent_child_table(parent_child_list)
newick_str = tree.write(format=1)
with open("tetra.tree", "w") as fh:
    fh.write(newick_str)

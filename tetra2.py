from ete3 import Tree
import numpy as np
import pandas as pd

class PTree:
    def __init__(self, df):
        self.df = df
        self.check_df()
        self.calc_branch_lengths()
        self.newick_string = self.generate_newick_string()

    def get_row(self, node):
        return self.df.loc[df['node'] == node, :]

    def get_parent(self, node):
        row = self.get_row(node)
        parent = row['parent'].iloc[0]
        return self.df.loc[self.df['node'] == parent, :]

    def get_children(self, node):
        return self.df.loc[self.df['parent'] == node, :]

    def _get_col(self, node, col):
        row = self.get_row(node)
        return row[col].iloc[0]

    def get_origin(self, node):
        return self._get_col(node, 'origin')

    def get_tip(self, node):
        return self._get_col(node, 'tip')

    def get_branch_length(self, node):
        return self._get_col(node, 'branch.length')

    def get_reference(self, node):
        return self._get_col(node, 'reference')

    def is_root(self, node):
        row = self.get_row(node)
        return row['parent'].iloc[0] is None

    def is_tip(self, node):
        row = self.get_row(node)
        return ~np.isnan(row['tip'].iloc[0])

    def has_parent(self, node):
        parent = self.get_parent(node)
        return parent.shape[0] > 0

    def has_children(self, node):
        children = self.get_children(node)
        return children.shape[0] > 0

    def check_df(self):

        nodes = self.df['node']
        if len(nodes) != len(nodes.drop_duplicates()):
            raise Exception('node names contain duplicates')
        if len(nodes) != len(nodes.dropna()):
            raise Exception('some node names are missing')

        parents = self.df['parent']
        if  len(parents.dropna()) == len(parents):
            raise Exception('no root node')
        if  len(parents.dropna()) < (len(parents) - 1):
            raise Exception(f'multiple root nodes:\n{self.df.loc[self.df['parent'].isna(), :]}')

        for node in nodes:  # node = 'Amniota'
            row = self.get_row(node)
            if self.is_root(node):
                if self.is_tip(node):
                    raise Exception(f'root should not be a tip:\n{row}')
                continue
            if not self.has_parent(node):
                raise Exception(f'parent of {node} does not exist:\n{row}')
            parent = self.get_parent(node)
            if self.get_origin(node) >= parent['origin'].iloc[0]:
                raise Exception(f'{node} should be younger than its parent:\n{pd.concat([row, parent])}')
            if self.has_children(node):
                if self.is_tip(node):
                    raise Exception(f'{node} is a tip but has children:\n{row}')
                children = self.get_children(node)
                if len(children['origin'].unique()) > 1:
                    raise Exception(f'{node} children have different origins:\n{children}')
                if self.get_origin(node) <= children['origin'].iloc[0]:
                    raise Exception(f'{node} should be older than its children:\n{pd.concat([row, children])}')

    def calc_branch_lengths(self):
        for node in df['node']:
            t0 = self.get_origin(node)
            if self.is_tip(node):
                t1 = self.get_tip(node)
            else:
                children = self.get_children(node)
                t1 = children['origin'].iloc[0]
            self.df.loc[self.df['node'] == node, 'branch.length'] = abs(t1 - t0)
            self.df = self.df[['node', 'parent', 'origin', 'tip', 'branch.length', 'reference']]

    def generate_newick_string(self):
        parent_child_list = self.df[['parent', 'node', 'branch.length']].apply(tuple, axis=1).tolist()
        tree = Tree.from_parent_child_table(parent_child_list)
        return tree.write(format=1)

    def write_newick_string(self, file):
        with open(file, "w") as fh:
            print(f'writing newick string to {file}')
            fh.write(self.newick_string)

    def __repr__(self):
        return self.df.to_string()

    def __str__(self):
        return self.df.to_string()


df = pd.DataFrame([
    # 	Carboniferous 358.9 - 298.9
    ('Amniota', None,  323.2, None, ''),
    ('Synapsida', 'Amniota', 320, None, ''),
    ('Sauropsida', 'Amniota', 320, None, ''),
    ('Diapsida', 'Sauropsida', 300, None, ''),
    ('Anapsida', 'Sauropsida', 300, 210, ''),
    # 	Permian 298.9 - 252
    # ('Testudines', 'Diapsida', 260, 0, ''),
    # 	Triassic 252 - 201.4
    ('Archosauria', 'Diapsida', 250, None, ''),
    ('Lepidosauria', 'Diapsida', 250, None, ''),
    # ('Thalattosauria', 'Diapsida', 300, 66, ''),
    # ('Ichthyosauromorpha', 'Diapsida', 300, 66, ''),
    # ('Sauropterygia', 'Sauropsida', 250, 66, ''),
    ('Therapsida', 'Synapsida', 279.5, None, ''),
    ('Sphenacodontids', 'Synapsida', 279.5, 252, ''),
    ('Cynodonts', 'Therapsida', 260, None, ''),
    ('Dicynodonts', 'Therapsida', 260, 201.4, ''),
    ('Avemetatarsalia', 'Archosauria', 245, None, ''),
    ('Pseudosuchia', 'Archosauria', 245, None, ''),
    # ('Placodonts', 'Sauropterygia', 250, 201.4, ''),
    # ('Nothosaurs', 'Sauropterygia', 250, 201.4, ''),
    # ('Plesiosaurs', 'Sauropterygia', 250, 66, ''),
    ('Dinosauria', 'Avemetatarsalia', 240, None, ''),
    ('Pterosauria', 'Avemetatarsalia', 240, 66, ''),
    ('Saurischia', 'Dinosauria', 235, None, ''),
    ('Ornithischia', 'Dinosauria', 235, 66, ''),
    ('Crocodylomorpha', 'Pseudosuchia', 225, None, ''),
    ('Aetosauria', 'Pseudosuchia', 225, 201.4, ''),
    ('Ornithosuchia', 'Pseudosuchia', 225, 201.4, ''),
    ('Rauisuchia', 'Pseudosuchia', 225, 201.4, ''),
    ('Mammalia', 'Cynodonts', 225, None, ''),
    # 	Jurassic 201.4 - 145
    ('Crocodilians', 'Crocodylomorpha', 200, 0, ''),
    # ('extinct-croc-relatives', 'Crocodylomorpha', 200, 66, ''),
    ('Monotremata', 'Mammalia', 161.5, 0, ''),
    ('Theria', 'Mammalia', 161.5, None, ''),
    ('Theropoda', 'Saurischia', 230, None, ''),
    ('Sauropodomorpha', 'Saurischia', 230, 66, ''),
    # ('ceratopsians', 'Ornithischia', 150, 66, ''),
    # ('pachycephalosaurs', 'Ornithischia', 150, 66, ''),
    # ('stegosaurs', 'Ornithischia', 150, 66, ''),
    # ('ankylosaurs', 'Ornithischia', 150, 66, ''),
    # 	Cretaceous 145 - 66
    ('Rhynchocephalia', 'Lepidosauria', 242, 0, ''),
    ('Squamata', 'Lepidosauria', 242, 0, ''),
    # ('Snakes', 'Lepidosauria', 145, 0, ''),
    # ('Lizards', 'Lepidosauria', 145, 0, ''),
    ('Metatheria', 'Theria', 110, 0, ''),
    ('Eutheria', 'Theria', 110, 0, ''),
    # ('Xenarthra', 'Eutheria', 58, 0, ''),
    # ('Afrotheria', 'Eutheria', 100, 0, ''),
    # ('Euarchontoglires', 'Eutheria', 100, 0, ''),
    # ('Laurasiatheria', 'Eutheria', 100, 0, ''),
    ('extinct-theropods', 'Theropoda', 150, 66, ''),
    ('birds', 'Theropoda', 150, 0, ''),
    # 	Paleogene 66 - 23.03
    # 	Neogene 23.03 - 2.58
    # 	Quaternary 2.58 - 0
], columns=['node', 'parent', 'origin', 'tip', 'reference'])
# lepidosaurs (lizards and snakes)
# print(df)

ptree = PTree(df)
print(ptree)
# print(ptree.newick_string)
ptree.write_newick_string('tetra.tree')

# node_id left_child_id right_child_id radius center_coordinates

import sys
import os


class Node:
    def __init__(self, node_id, left_child_id, right_child_id, radius, center, tree):
        self.node_id = node_id
        self.left_child_id = left_child_id
        self.right_child_id = right_child_id
        self.radius = radius
        self.center = center
        self.tree = tree

    def __eq__(self, other):
        result = None

        if self.radius != other.radius or self.center != other.center:
            result = False
        elif self.left_child_id == -1 and self.left_child_id != other.left_child_id:
            result = False
        elif self.right_child_id == -1 and self.right_child_id != other.right_child_id:
            result = False

        elif self.right_child_id == -1 and self.right_child_id == other.right_child_id and self.left_child_id == -1 and self.left_child_id == other.left_child_id:
            result = True

        if result == False:
            print(self.node_id, self.left_child_id,
                  self.right_child_id, self.radius, self.center)
            print(other.node_id, other.left_child_id,
                  other.right_child_id, other.radius, other.center)
        if not result is None:
            return result

        return (self.left_child_id == -1 or other.left_child_id == -1 or self.tree.nodes[self.left_child_id] == other.tree.nodes[other.left_child_id]) and\
            (self.right_child_id == -1 or other.right_child_id == -
             1 or self.tree.nodes[self.right_child_id] == other.tree.nodes[other.right_child_id])


class Tree:
    def __init__(self, num_dims, num_nodes):
        self.num_dims = num_dims
        self.num_nodes = num_nodes
        self.nodes = {}

    def __eq__(self, other):
        return self.nodes[0] == other.nodes[0]


def readFile(file):
    with open(file) as f:
        specs = f.readline()
        num_dims, num_nodes = specs.split(" ")
        tree = Tree(num_dims, num_nodes)
        while f.readable():
            line = f.readline().strip()
            if line is None or line == "":
                break
            line = line.split(" ")
            node_id = int(line[0])
            left_child_id = int(line[1])
            right_child_id = int(line[2])
            radius = line[3]
            center = " ".join(line[4:])
            tree.nodes[node_id] = Node(
                node_id, left_child_id, right_child_id, radius, center, tree)

    return tree


def main():
    if len(sys.argv) != 3:
        print("Usage: python verifier.py <file1.tree> <file2.tree>")
        return -1

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    if not os.path.exists(file1):
        print("File {} not found".format(file1))
        return -1

    if not os.path.exists(file2):
        print("File {} not found".format(file2))
        return -1

    tree = readFile(file1)
    tree2 = readFile(file2)
    print("Are they equivalent?", tree == tree2)


if __name__ == "__main__":
    main()

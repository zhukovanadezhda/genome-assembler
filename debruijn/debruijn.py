#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
from pathlib import Path
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
    random_layout,
    draw,
    spring_layout,
    isolates,
    shortest_path
)
import matplotlib
from operator import itemgetter
import random

random.seed(9001)
import itertools
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
from typing import Iterator, Dict, List

matplotlib.use("Agg")

__author__ = "Nadezhda Zhukova"
__copyright__ = "Université Paris Cité"
__credits__ = ["Nadezhda Zhukova"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Nadezhda Zhukova"
__email__ = "nadiajuckova@gmail.com"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i", 
        dest="fastq_file", 
        type=isfile, 
        required=True, 
        help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", 
        type=int, 
        default=22, 
        help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", 
        dest="graphimg_file", 
        type=Path, 
        help="Save graph as an image (png)"
    )
    return parser.parse_args()


def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterates over the read sequences.
    """
    with fastq_file.open('r') as file:
        for header in file:
            sequence = next(file).strip()  # Read sequence
            next(file)  # Skip plus line
            next(file)  # Skip quality line
            yield sequence


def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers of size kmer_size.
    """
    for i in range(0, len(read) - kmer_size + 1):
        yield read[i : i + kmer_size]


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph of all kmer substring and weight (occurrence).
    """
    graph = DiGraph()
    for kmer, weight in kmer_dict.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=weight)
    return graph


def remove_paths(graph: DiGraph, path_list: List[List[str]],
                 delete_entry_node: bool, delete_sink_node: bool) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected nodes in
    the graph.

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of paths
    :param delete_entry_node: (boolean) True -> Remove the first node of a path
    :param delete_sink_node: (boolean) True -> Remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        for i in range(len(path) - 1):
            if delete_entry_node and i == 0:
                graph.remove_node(path[i])
            elif delete_sink_node and i == len(path) - 2:
                graph.remove_node(path[i + 1])
            else:
                graph.remove_edge(path[i], path[i + 1])

    # Remove isolated nodes
    isolated_nodes = list(isolates(graph))
    graph.remove_nodes_from(isolated_nodes)

    return graph


def select_best_path(graph: DiGraph, path_list: List[List[str]],
                     path_length: List[int], weight_avg_list: List[float],
                     delete_entry_node: bool = False,
                     delete_sink_node: bool = False) -> DiGraph:
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of paths
    :param path_length: (list) A list of lengths of each path
    :param weight_avg_list: (list) A list of average weights of each path
    :param delete_entry_node: (bool) True -> Remove the first node of a path
    :param delete_sink_node: (bool) True -> Remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    # Compare based on stdev of the average weight
    if len(weight_avg_list) > 1 and statistics.stdev(weight_avg_list) > 0:
        # Select path with the highest average weight
        best_path_index = weight_avg_list.index(max(weight_avg_list))
    elif len(path_length) > 1 and statistics.stdev(path_length) > 0:
        # If weights are equal, compare path lengths
        best_path_index = path_length.index(max(path_length))
    else:
        # If both weights and lengths are equal, choose randomly
        best_path_index = randint(0, len(path_list) - 1)

    # Remove non-best paths
    paths_to_remove = [path for i, path in enumerate(path_list) 
                       if i != best_path_index]
    
    graph = remove_paths(graph, paths_to_remove, 
                         delete_entry_node, 
                         delete_sink_node)

    return graph


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, 
                 ancestor_node: str, 
                 descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue by selecting the best path.

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # Find all simple paths between the ancestor and descendant
    all_paths = list(all_simple_paths(graph, 
                                      source=ancestor_node, 
                                      target=descendant_node))
    # Compute the average weight of the path
    weight_avgs = [path_average_weight(graph, path) for path in all_paths]
    # Compute the length of each path
    path_lengths = [len(path) for path in all_paths]
    # Choose the best path and remove others
    graph = select_best_path(graph, all_paths, path_lengths, weight_avgs)
    
    return graph


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False 
    
    # Iterate through each node to check for bubbles
    for node in graph.nodes():
        # Get the list of predecessors for the current node
        liste_predecesseurs = list(graph.predecessors(node))
        
        # If there is more than one predecessor, a potential bubble exists
        if len(liste_predecesseurs) > 1:
            # Check each unique combination of predecessors (i, j)
            for i, j in itertools.combinations(liste_predecesseurs, 2):
                # Find the lowest common ancestor of the two predecessors
                noeud_ancetre = lowest_common_ancestor(graph, i, j)
                # If a common ancestor is found, a bubble exists
                if noeud_ancetre is not None:
                    bubble = True
                    break
            if bubble:
                break
    
    # Recursive simplification of the graph
    if bubble:
        # Solve the bubble between the ancestor and the current node
        graph = simplify_bubbles(solve_bubble(graph, noeud_ancetre, node))
    
    return graph


def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    pass


def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    pass


def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors
    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    return [node for node in graph.nodes() 
            if not list(graph.predecessors(node))]


def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    return [node for node in graph.nodes() 
            if not list(graph.successors(node))]


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            # Get all simple paths from the start node to the end node
            paths = all_simple_paths(graph, start_node, end_node)
            for path in paths:
                # Build the contig from the path
                # Start with the first k-mer
                contig = path[0]
                for node in path[1:]:
                    # Append only the last nucleotide of each subsequent k-mer
                    contig += node[-1]  
                # Store contig and its length
                contigs.append((contig, len(contig)))  
    
    return contigs


def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with output_file.open("w", newline='\n') as file:
        for i, (contig, length) in enumerate(contigs_list):
            file.write(f">contig_{i} len={length}\n")
            file.write(textwrap.fill(contig, 80) + "\n")


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None: # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) 
              in graph.edges(data=True) if d["weight"] > 3]
    # print(elarge)
    esmall = [(u, v) for (u, v, d) 
              in graph.edges(data=True) if d["weight"] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, 
        edge_color="b", style="dashed"
    )
    # nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file.resolve())


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == "__main__":  # pragma: no cover
    main()

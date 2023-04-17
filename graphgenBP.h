#ifndef _graphgenBP_h
#define _graphgenBP_h



typedef struct /* the bipartite graph data structure */
{
	long n; // numver of vertices in both sides
	long nrows; // number of vertices in the left side
	long m; // number of edges
	long* vtx_pointer; // an array of size n+1 storing the pointer in endV array
    long* endV; //an array of size m that stores the second vertex of an edge.
	double* weight; // not used in unweighted graph
} graph;



void process_mtx_compressed(char *fname, graph* bGraph);
void fast_mtx_read_build(char *fname, graph* bGraph);
bool isEqual(graph* bGraph1, graph* bGraph2);

/** Clean up. */
void free_graph (graph* bGraph);
graph* swap_side(graph* bGraph);

#endif

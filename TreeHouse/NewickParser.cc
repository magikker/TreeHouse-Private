//file: NewickParser.cc
#include "NewickParser.h"
#include <../gulo/Phylo/ExtendedNewickFormat>
#include <../gulo/Graph/graphparser.hpp>

using namespace Gulo;

/*
class custom_dfs_visitor : public boost::default_dfs_visitor { 
	public: template < typename Vertex, typename Graph >
	void discover_vertex(Vertex u, const Graph & g)	const { std::cout << "At " << u << std::endl; }
	template < typename Edge, typename Graph >
	void examine_edge(Edge e, const Graph& g) const { std::cout << "Examining edges " << e << std::endl;} 
}; 
*/

struct BGraph
{  
    typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
    typedef boost::property<boost::vertex_name_t, std::string> VertexNameProperty;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexNameProperty, EdgeWeightProperty> AdjacencyList;
    typedef boost::graph_traits<AdjacencyList>::vertex_descriptor Node;
    typedef boost::graph_traits<AdjacencyList>::edge_descriptor Edge;

    BGraph():
    edge_length(get(boost::edge_weight, g)),
    node_label(get(boost::vertex_name, g))
    {
    }

    AdjacencyList g;
    boost::property_map<AdjacencyList, boost::edge_weight_t>::type edge_length;
    boost::property_map<AdjacencyList, boost::vertex_name_t>::type node_label;
    
    Node root;
};

/*
struct BGraph
{  
    typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
    //typedef boost::property<boost::vertex_name_t, std::string> VertexNameProperty;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, node_name, EdgeWeightProperty> AdjacencyList;
    typedef boost::graph_traits<AdjacencyList>::vertex_descriptor Node;
    typedef boost::graph_traits<AdjacencyList>::edge_descriptor Edge;

    BGraph():
    edge_length(get(boost::edge_weight, g)),
    node_label(get(&node_name::name, g))
    {
    }

    AdjacencyList g;
    boost::property_map<AdjacencyList, boost::edge_weight_t>::type edge_length;
    boost::property_map<AdjacencyList, get(&node_name::name, g)>::type node_label;
    
    Node root;
};
*/
/*
struct children_counting_dfs_visitor : public boost::default_dfs_visitor { 

    clade_collector_dfs_visitor(BGraph * g)
    : graph(g)
    {
		for(size_t i = 0; i < boost::num_vertices(graph->g); i++){
			num_children.push_back(0);
		}
    }

	
	//public: template < typename Vertex, typename Graph >
	//void discover_vertex(Vertex u, const Graph & g)	const { 
	//	std::cout << "At " << u << std::endl;
	//}
	
	
	public: template < typename Vertex, typename Graph >
	void finish_vertex(Vertex u, const Graph & g) { 
		cout << graph->node_label[u] << endl; 
		typedef boost::graph_traits<Graph> GraphTraits;

		typename boost::property_map<Graph, boost::vertex_index_t>::type index = get(boost::vertex_index, g);

		if (out_degree(u,g) == 0){
			//leaf node. 
			num_children[index[u]] = 0;
		}
		else{
			typename GraphTraits::out_edge_iterator out_i, out_end;
			typename GraphTraits::edge_descriptor e;
			//add the count for each child to the current node
			int grandkids = 0;
			for (tie(out_i, out_end) = out_edges(u, g); out_i != out_end; ++out_i) {
				e = *out_i;
				Vertex targ = target(e, g);
				grandkids += num_children[index[targ]];
			}
			num_children[index[u]] = grandkids + out_degree(u,g);
			
		}
		std::cout << "At " << u << " there are " << num_children[index[u]] << " children." << std::endl;
	}
	

    BGraph * graph;
    vector<int> num_children;
};  
*/


struct clade_collector_dfs_visitor : public boost::default_dfs_visitor { 

    clade_collector_dfs_visitor(BGraph * g, vector<set <unsigned int>> & tc)
    : graph(g), treeclades(tc){
    }

	public: template < typename Vertex, typename Graph >
	void finish_vertex(Vertex u, const Graph & g) { 
		//cout << graph->node_label[u] << endl; 
		typedef boost::graph_traits<Graph> GraphTraits;

		typename boost::property_map<Graph, boost::vertex_index_t>::type index = get(boost::vertex_index, g);

		//leaf node // trivial clade. 
		if (out_degree(u,g) == 0){
			//add to label map and get position 
			treeclades[index[u]].insert(::biparttable.lm.add(graph->node_label[u]));
		}
		else{
			typename GraphTraits::out_edge_iterator out_i, out_end;
			typename GraphTraits::edge_descriptor e;
			//add the clade members of each child to the current node's clade
			for (tie(out_i, out_end) = out_edges(u, g); out_i != out_end; ++out_i) {
				e = *out_i;
				Vertex targ = target(e, g);
				treeclades[index[u]].insert(treeclades[index[targ]].begin(), treeclades[index[targ]].end() );
			}
		}
		//cout << "treeclades[index[u]].size() = "<< treeclades[index[u]].size() << endl;
		//for(size_t i = 0; i < treeclades[index[u]].size(); i++){
		//		cout << treeclades[index[u]][i] << " ";
		//}
		//cout << endl;		
	}
	
	/*
	template < typename Edge, typename Graph >
	void examine_edge(Edge e, const Graph& g) const { std::cout << "Examining edges " << e << std::endl;}
	*/
    BGraph * graph;
    vector<set <unsigned int>> & treeclades;
    //vector<vector <unsigned int>> treeclades;
};  



struct MyVisitor : public boost::default_dfs_visitor {
    MyVisitor(BGraph * g)
    : graph(g)
    {
    }

    template <class Edge, class Graph>
    void examine_edge(Edge e, const Graph & g)
    {
        std::cout << graph->node_label[source(e, g)]
                << "->"
                << graph->node_label[target(e,g)] << std::endl;
    }

    BGraph * graph;
};

namespace Gulo
{
    template <>
    struct GraphParserMapping<ExtendedNewickFormat, BGraph>
    {
        typedef BGraph::Node VertexIdentifier;
        typedef BGraph::Edge EdgeIdentifier;

        VertexIdentifier addVertex(BGraph & graph, const std::vector<std::string> & comments)
        {
            return add_vertex(graph.g);
        }

        EdgeIdentifier addEdge(BGraph & graph, VertexIdentifier from, VertexIdentifier to)
        {
            return (add_edge(from, to, graph.g)).first;
        }

        EdgeIdentifier trailingEdge(BGraph&)
        {
            return EdgeIdentifier();
        }

        void setLength(BGraph & graph, EdgeIdentifier edge, double length, const std::vector<std::string> & comments)
        {
            if (edge != EdgeIdentifier())
                graph.edge_length[edge] = length;
        }

        void setLabel(BGraph & graph, VertexIdentifier node, std::string && label, const std::vector<std::string> & comments)
        {
            graph.node_label[node] = std::move(label);

        }

        void setHybridVertexType(BGraph & g, VertexIdentifier n, const std::string & label, const std::vector<std::string> & comments)
        {
            // i don't really care about the hybrid node type  
        }

        void setRoot(BGraph & g, VertexIdentifier n)
        {
            g.root = n;
        }

        void setUnrooted(BGraph & g)
        {
            g.root = VertexIdentifier();
        }

        void abort(BGraph & g)
        {
            g.g.clear();
        }
    };
}

// This is a visitor to traverse the generated graph, it's boost code


struct less_than_key
{
    inline bool operator() (const set <unsigned int> & set1, const set <unsigned int> & set2)
    {
        return (set1.size() < set2.size());
    }
};




void build_tree( vector<set <unsigned int>> treeclades, BGraph graph ){
    typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
    typedef boost::property<boost::vertex_name_t, std::string> VertexNameProperty;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexNameProperty, EdgeWeightProperty> AdjacencyList;
    typedef boost::graph_traits<AdjacencyList>::vertex_descriptor Node;
    typedef boost::graph_traits<AdjacencyList>::edge_descriptor Edge;
	std::sort(treeclades.begin(), treeclades.end(), less_than_key());
	map<unsigned int, Node > temp;

	for(vector<set<unsigned int>>::iterator it = treeclades.begin(); it != treeclades.end(); ++it){
		if((*it).size() == 1){
			for(auto f : *it) {
				temp[f] = boost::add_vertex(graph.g);
			} 
		}
	}


}




int graphtest() {

    //std::string newick = "(clouded_leopard,(snow_leopard,(tiger,(jaguar,lion,leopard))));";

    #include <fstream>
    //string fileName = "test/data/panthera.tre";
	string fileName = "../../synthesis_trees/life_synth.tre";
	ifstream inFile;
	inFile.open(fileName);
	if (inFile.fail()) {
      cerr << "unable to open file "<< fileName <<" for reading" << endl;
      exit(1);
	}

	map<set<unsigned int>, vector<unsigned int>> tempstructure;
	map<set<unsigned int>, vector<unsigned int>>::iterator iter;
	std::string line;
	unsigned int count = ::biparttable.NumTrees + 1;
	while (std::getline(inFile, line)){
	    ::biparttable.NumTrees +=1;
		BGraph graph;

		GraphParser<ExtendedNewickFormat, BGraph> parser;
		auto result = parser.parse(line, graph);
		if (result.hasError()) {
			std::cout << result.errorMessage() << std::endl << result.errorSnippet() << std::endl;
			return -1;
		}
		vector<set <unsigned int>> treeclades;
		treeclades.resize(num_vertices(graph.g));
		clade_collector_dfs_visitor visc(&graph, treeclades);
		
	    //std::cout << "The graph contains " << num_vertices(graph.g) << " nodes" << std::endl;
	    
	    boost::depth_first_search(graph.g, visitor(visc));
	    //std::cout << "The graph contains " << num_vertices(graph.g) << " nodes" << std::endl;
	    
	    //for(size_t i = 0; i < treeclades.size(); i++){
		//	for(size_t j = 0; j < treeclades[i].size(); j++){
		//		cout << treeclades[i][j] << " ";
		//	}
		//	cout << endl;
		//}
		for(unsigned int i = 0; i < treeclades.size(); i++){
			
			iter = tempstructure.find(treeclades[i]);
			if (iter != tempstructure.end() ){
				iter->second.push_back(count);
			}
			else{
				tempstructure[treeclades[i]] = vector<unsigned int>({count});
			}	
		}
		count++;
	}
	
    for (std::map<set<unsigned int>, vector<unsigned int>>::iterator it=tempstructure.begin(); it!=tempstructure.end(); ++it){
		std::cout << it->first.size() << " => " << it->second.size() << '\n';
	}

	cout << tempstructure.size() << endl;
	cout << sizeof(tempstructure) << endl;
	cout << sizeof(::biparttable.CladeMap) << endl;
    //auto result = parser.parse(newick, graph);
    //if (result.hasError()) {
    //    std::cout << result.errorMessage() << std::endl << result.errorSnippet() << std::endl;
    //    return -1;
    //}

    ////////std::cout << "The graph contains " << num_vertices(graph.g) << " nodes" << std::endl;

    //MyVisitor vis(&graph);
    //std::vector < boost::default_color_type > color_map(boost::num_vertices(graph.g));
    //depth_first_visit(graph.g,
    //                graph.root,
    //                vis,
    //                make_iterator_property_map(color_map.begin(), get(boost::vertex_index, graph.g), color_map[0]));

    ///////////////clade_collector_dfs_visitor visc(&graph);
    //depth_first_visit(graph.g,
    //                graph.root,
    //                visc,
    //                make_iterator_property_map(color_map.begin(), get(boost::vertex_index, graph.g), color_map[0]));

	//////////boost::depth_first_search(graph.g, visitor(visc));


/*
 typedef boost::property<boost::edge_weight_t, int>  EdgeWeightProperty; 
 typedef boost::adjacency_list  < boost::listS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty>  mygraph;

 mygraph g; 
 boost::add_edge (0, 1, 8, g);
 boost::add_edge (0, 3, 18, g);
 boost::add_edge (1, 2, 20, g);
 boost::add_edge (2, 3, 2, g);
 boost::add_edge (3, 1, 1, g);
 boost::add_edge (1, 3, 7, g);
 custom_dfs_visitor vis;
 boost::depth_first_search(g, visitor(vis));
*/

 return 1;
}


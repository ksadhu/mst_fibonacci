/*
   To run the program:
compile:
	g++ mst.cpp
	----> this generated a.out in current directory
run:
	./a.out -s input 
	-----> for simple scheme with input as input file
	./a.out -f input 
	-----> for fibonacci scheme with input as input file
	./a.out -r n d 
	-----> for random mode n vertices and d% density.

   This program implements prims algorithm for minumim spannig trees in two ways.
   1.Simple sheme: In this scheme the weights extracted and are updated in an array
   2.Fibonacci heap: In this scheme the weights are extracted and are updated through a fibonacci heap
 	The code supports two modes for running the program.
	1.Random mode
	2.User input mode
	In random mode, a graph is generated with number of vertices and density as input. Using rand fuction, all the edges and their
	costs are generated and added to the graph. Then both the schemes are called to run prims algorithm.
	In user input mode, the input from the a file is read and along with the input as which scheme should the algorithm run in.
	Author: Kumar Sadhu 
	UFID:29493662
   This code has been implemented by me for most of the part except for a small portion of the The graph representation of through 
   adjacency list which has been taken from a reference on geeksforgeeks as a reuse package in the code. 
   The implementation of Fibonacci heap is based on the algorithms described on CLRS.
*/

#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAX_WEIGHT 2000
using namespace std;
//This is the node structure for a node in fibonacci heap
typedef struct f_heap_node{
	int degree;       // Number of children of this node
	bool c_cut; // Whether this node is marked
	struct f_heap_node* next;   // Next and previous elements in the list
	struct f_heap_node* prev;   
	struct f_heap_node* parent; // Parent of this node, if any.
	struct f_heap_node* child;  // Child node, if any.
	int vertex;     // Element being stored here
	int weight; // Its priority
}fnode;
// This is the fibonacci heap data structure which has min pointer
// of type fibonacci node and total nodes in the heap
typedef struct f_heap{
	fnode* minNode;
	int nodes;
}fheap;
// This function creates an empty heap and initializes its elements
// It is used only when we are creating a new f-heap
fheap* makeheap(fheap* &root){
	root= (fheap*)malloc(sizeof(fheap));
	root->minNode=NULL;
	root->nodes=0;
	return root;
}
// This function combines the circular doubly linked lists
// n1 and n2 and returns a pointer to the merged doubly linked list.
fnode* mergeroots(fnode* &n1,fnode* &n2){
	// base case 1, when both are null return null.
	if(n1==NULL and n2==NULL)
		return NULL;
	// base case 2, when either of them are null return the other.
	else if(n1==NULL and n2!=NULL)
		return n2;
	else if(n2==NULL and n1!=NULL)
		return n1;
	// case 3: when both are not null, we need to manipulate pointers for
	// n1 and n2 such that we get a single pointer for entire list.
	else{
		// declare a temporary pointer to store next of n1
		fnode* temp = n1->next;
		// assign next of n2 to next of n1. So we are breaking n1->next
		n1->next = n2->next;
		// this assigns n1 as the prev pointer to above node. we are half done.
		n1->next->prev = n1;
		// assign the original lost pointer to the next of n2
		n2->next = temp;
		// this step connects n2 to the n1 back. 
		n2->next->prev = n2;
		// here we return the pointer which ever is minimum of n1 and n2
		if(n1->weight<n2->weight)
			return n1;
		else
			return n2;
	}
}
// This method inserts a fibonacci node into fibonacci heap h with values vertex and weight
fnode* finsert(fheap* &h,fnode* &node,int vertex,int weight){
	// node creation and initialization of node fields.
	node= (fnode*)malloc(sizeof(fnode));
	node->degree=0;
	node->parent=NULL;
	node->child=NULL;
	node->next=node;
	node->prev=node;
	node->c_cut=false;
	node->vertex=vertex;
	node->weight=weight;
	// get current minpointer in heap
	fnode* cur_min=h->minNode;
	// mergewith newly created node
	cur_min=mergeroots(cur_min,node);
	// we again have the min pointer set for our heap
	h->minNode=cur_min;
	//h->minNode=mergeroots(h->minNode,node);
	//if(h->minNode==NULL or ((h->minNode)->weight)>node->weight)
	//	h->minNode=node;
	// increase the count of nodes in the heap
	(h->nodes)++;
	// return a pointer to the newly created node.
	// we make use of this pointer while running prims algorithm to locate a node
	// in a heap.
	return node;
}
// This method merges two fibonacci heaps h1 and h2 into one and returns merged one.
// although this method is not called anywhere in the code, its good to have it implemented!!
fheap* mergeheaps(fheap* h1,fheap* h2){
	fheap* ans;
	// create a new heap
	ans=makeheap(ans);
	ans->minNode=h1->minNode;
	// merge minnodes of both h1 and h2
	ans->minNode=mergeroots(ans->minNode,h2->minNode);
	// assign overall minnode
	if(h1->minNode==NULL or (h2->minNode!=NULL and h2->minNode<h1->minNode))
		ans->minNode=h2->minNode;
	// set nodes in new heap to sum of nodes in h1 and h2
	ans->nodes=h1->nodes+h2->nodes;
	// destroy old heaps h1 and h2
	free(h1);
	free(h2);
	// return new heap
	return ans;
}
// This method removes the node chld from its parent prnt in the heap h 
void cut(fheap* &h ,fnode* &chld ,fnode* &prnt ){
	//if chld has siblings,then remove it from sibling list by traversing
	if(chld->next!=chld){
		chld->next->prev=chld->prev;
		chld->prev->next=chld->next;
	}
	//if the node being removed is the child node of parent, we have to set
	// some other sibling of chld as child for prnt. 
	if(prnt->child == chld){
		if (chld->next != chld) {
			prnt->child = chld->next;
		}
		// if its the only child then prnt has child nomore so set to null.
		else
			prnt->child=NULL;
	}
	//removing the child means we have to decrease the degree of parent
	(prnt->degree)--;
	// make the child as a singleton node
	chld->next=chld->prev=chld;
	// merge it with rootlists i.e with min pointer node
	h->minNode=mergeroots(h->minNode,chld);
	// set its parent to null and child cut value to false.
	chld->parent=NULL;
	chld->c_cut=false;
}
// This recursive function removes the nodes from heap if the child cut
// values are true. It goes from child to parent until it sees a parent
// whose child cut value is false. It sets the child cut of last parent 
// as true since it just lost a child.
void cascade(fheap* &h,fnode* &y){
	fnode* z=y->parent;
	if(z!=NULL){
		// case when parent is having child cut false
		if(y->c_cut==false){
			y->c_cut=true;
		}
		// remove this node and recurse up the path until child cut is false
		else{
			cut(h,y,z);
			cascade(h,z);
		}
	}
}
// This method sets the weight field of node n to amount in the heap h
// Based on the new value it may be removed from the parent and called for
// cut and cascade methods.
void decrease_weight(fheap* &h,fnode* &n,int amount){
	if(amount>n->weight){
		cout<<"Error! New key is greater than current"<<endl;
	}
	else{
		// set the new weight for the node n
		n->weight=amount;
		fnode* prnt=n->parent;
		// new weight is less than parent weight, so it must be cut 
		// from the parent
		if(prnt!=NULL and n->weight<prnt->weight){
			// we have to cut it!
			cut(h,n,prnt);
			cascade(h,prnt);
		}
		//just update min pointer
		if(n->weight<(h->minNode)->weight)
			h->minNode=n;
	}
}
// This method is the crucial method for fibonacci heaps where all the actual time
// is spent in pairwise merging of all the trees, to ensure at the end of the operation
// there are no two trees with the same degree.
// we scan the root list with the help of min pointer in the heap and combine trees with
// same degree into a single tree. We do this with the help of a auxillary table to see if
// there is a tree with degree already existing in our heap.
void consolidate(fheap* &h){
	// trees to scan in the heap.
	vector<fnode*> heaptrees;
	fnode* temp=h->minNode;
	fnode* aux=temp;
	if(temp!=NULL){
		// push all the trees to the table so we can scan through them later.
		while(heaptrees.empty() or temp!=aux){
			if(temp!=NULL){
				heaptrees.push_back(temp);
				temp=temp->next;
			}
		}
		// Auxillary table to store node address with each degree. While scanning
		// through the heaptrees we see if there exists a tree with degree same as current
		// tree. If yes, they are merged and table is updated with new degrees.
		vector<fnode*> table;
		for(vector<fnode*>::iterator it = heaptrees.begin() ; it != heaptrees.end(); ++it){
			fnode* temp=*it;
			while(true){
				while (temp->degree >= table.size()){
					table.push_back(NULL);
				}
				// if there is no tree waiting, break the loop
				if ((table[temp->degree]) == NULL) {
					table[temp->degree]=temp;
					break;
				}
				// else we have to pair wise merge.
				fnode* mrg=table[temp->degree];
				table[temp->degree]=NULL;
				// assign small and large nodes based on comparison of weights
				fnode* small = (mrg->weight<temp->weight)?mrg:temp;
				fnode* large = (mrg->weight<temp->weight)?temp:mrg;
				large->next->prev=large->prev;
				large->prev->next=large->next;
				large->next=large->prev=large;
				// merge larger one with smaller
				small->child=mergeroots(small->child,large);
				large->parent=small;
				// set the child cut as false for large node as its just added as a child.
				large->c_cut=false;
				// increment the degree of small node
				(small->degree)++;
				temp=small;
			}
			if(temp->weight<=(h->minNode)->weight){
				h->minNode=temp;
			}
		}
		//heaptrees.erase(heaptrees.begin(),heaptrees.end());
		//table.erase(table.begin(),table.end());
	}
}
// This method removes the minNode from the heap h and returns it. This is the most
// used function in prims algorithm as we need edge with minimum weight in each iteration.
fnode* removeMin(fheap* &h){
	// get the min node from heap
	fnode* ans=h->minNode;
	// decrease the node number from heap
	(h->nodes)--;
	// remove this node from the root list
	// if its a single node in the root list, we have our heap empty so set it to null
	if ((h->minNode)->next == h->minNode) {
		h->minNode = NULL;
	}
	// if there are multiple nodes in the root level, delete current node by manipulating 
	// pointers in the doubly linked list. Assign minNode for heap arbitrarily.
	else{
		(h->minNode)->prev->next = (h->minNode)->next;
		(h->minNode)->next->prev = (h->minNode)->prev;
		h->minNode = (h->minNode)->next;
	}
	// For all the children we need to set the parent as null as they will be promoted to
	// rootlevel list.
	if(ans->child!=NULL){
		fnode* chld=ans->child;
		do{
			chld->parent=NULL;
			chld=chld->next;
		}
		while(chld!=ans->child);
	}
	// merge the child list of deleted node with min pointer in the heap.
	h->minNode=mergeroots(h->minNode,ans->child);
	if(ans==NULL){
		return NULL;
	}
	// pairwise merge the trees
	consolidate(h);
	return ans;
}
// This is a representation of graph node.
struct AdjListNode{
	int dest;// The destination vertex j in the edge(i,j)
	int weight;// weight of edge (i,j)
	struct AdjListNode* next;// pointer to next node
};
// A struct to store the head pointer for each adjacency list
struct AdjList{
	struct AdjListNode *head;  // pointer to head node of list
};
// Graph contains total vertices V and an array of adjacency lists
struct Graph{
	int V;
	struct AdjList* array;
};
// constructs a new node with destination dest and next as null and returns it
struct AdjListNode* newAdjListNode(int dest){
	struct AdjListNode* newNode =
		(struct AdjListNode*) malloc(sizeof(struct AdjListNode));
	newNode->dest = dest;
	newNode->next = NULL;
	return newNode;
}
// This method allocates memory for the graph and initializes its fields
// and returns the created graph
struct Graph* createGraph(int V){
	struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
	graph->V = V;
	graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));
	int i;
	for (i = 0; i < V; ++i)
		graph->array[i].head = NULL;
	return graph;
}
// This method adds an edge with given parameters in the gragh strucute.
void addEdge(struct Graph* graph, int src, int dest, int weight){
	// create a dest node in the adjacency list
	struct AdjListNode* newNode = newAdjListNode(dest);	
	newNode->next = graph->array[src].head;
	newNode->weight=weight;
	//add it to the head pointer in the index of source vertex.
	graph->array[src].head = newNode;
	// as its undirected graph, it means we have to create another entry
	// in the dest adjacency list for the source node.
	newNode = newAdjListNode(src);
	newNode->next = graph->array[dest].head;
	newNode->weight=weight;
	graph->array[dest].head = newNode;
}
// Method which generated a random graph with input parameters vertices and
// density. Total edges are calculated.
struct Graph* generateRandomGraph(int vertices, int density){
	// calculate total edges based on n and d
	// generate triplet i j and k
	int num_nodes;
	num_nodes=vertices;
	// create and initialize a graph
	struct Graph* graph = createGraph(num_nodes);
	// an auxilary map to check if the edge generated is not duplicate
	// it is of the form (i,j) as a pair with 1 as the key.
	map< pair<int,int>,int> check_map;
	// calculate number of edges
	long double num_edges=num_nodes*(num_nodes-1)/2;
	num_edges=num_edges*density/100.0;
	cout<<"edges are : "<<num_edges<<endl;
	int count=0;
	// keep generating edges until the count is reached
	while(count<num_edges){
		int start,end,weight;
		start=rand()%num_nodes;
		end=rand()%num_nodes;
		// make sure always i is less than j in (i,j)
		if(start>end){
			int temp=start;
			start=end;
			end=temp;
		}
		// create edge pair and generate weight for edge
		pair<int, int> edge(start,end);
		weight=(rand()%1000)+1;
		// if its not yet inserted into graph insert it into graph as well as the map
		if(check_map.find(edge)==check_map.end() and start!=end){
			check_map.insert(make_pair(make_pair(start,end),1));
			//cout<<start<<" "<<end<<" "<<weight<<endl;
			addEdge(graph,start,end,weight);
			count++;
		}
	}
	return graph;
}
// A dfs routine to check if the graph is connected or not.
int dfs_count=1;
void dfs(struct Graph* g, int v, bool marked[]){
	marked[v]=true;	
	struct AdjListNode* tr_list = g->array[v].head;
	//printf("\n Adjacency list of vertex %d\n head ", v);
	while (tr_list){
	      if(marked[tr_list->dest]==false){
		      dfs_count++;
		      //if(dfs_count==g->V)
			//      return ;
		      dfs(g,tr_list->dest,marked);
	      }
	      //printf("-> %d(%d)", tr_list->dest,tr_list->weight);
	      tr_list = tr_list->next;
	}
}
// A helper method to check if graph is connected or not. It calls dfs internally
// and returns true if all vertices are visited in dfs and false otherwise.
bool checkConnected(struct Graph* g){
	int nodes=g->V;
	bool marked[nodes];
	for(int i=0;i<nodes;i++)
		marked[i]=false;
	dfs(g,0,marked);
	if(dfs_count==nodes){
		//cout<<"yay! connected!"<<endl;
		return true;
	}
	else{
		//cout<<"not connected!"<<endl;
	}
	// if there is any vertex which has not been visted after dfs return false.
	for(int i=0;i<nodes;i++){
		if(marked[i]==false)
			return false;
	}
}
// Prims implementaion using fibonacci scheme(as a priority queue to decide least weight).
clock_t prims_fscheme(struct Graph* g){
	//initialize fibonacci heap with vertices and weights as inf.
	clock_t Start, Time;
	// start the timer.
	Start = clock();
	//create heap and set nodes to number of vertices.
	fheap* root;
	root=makeheap(root);
	fnode* f1; 
	int nodes=g->V;
	int count=0,cost=0;
	// auxillary arrays to mark parents and check if the node is visited or not
	int parent[nodes];
	bool visited[nodes];
	list<int> output;
	// auxillary map to store the address of the node based on vertex number as key.
	map<int,fnode*> hmap;
	// insert all the vertices into f-heap with max weights and parent as -1
	// and mark visited as false.
	while(count<nodes){
		f1=finsert(root,f1,count,MAX_WEIGHT);
		hmap.insert(make_pair(count,f1));
		parent[count]=-1;
		visited[count]=false;
		count++;
	}
	int total=0;
	// set the weight of root to 0 to start the algorithm.
	root->minNode->weight=0;
	// while there are elements in the heap.
	while(root->nodes>0){
		// get a min element
		fnode* min=removeMin(root);
		// mark it as visited
		visited[min->vertex]=true;
		// traverse the adjacency list of this vertex.
		struct AdjListNode* tr_list = g->array[min->vertex].head;
		while (tr_list){
			// find the address of adjacent node in the heap
			fnode* out=hmap[tr_list->dest];
			//check if unvisited and its weight is more than edge weight
			if(visited[tr_list->dest]==false and out->weight>tr_list->weight){
				// update weight of this vertex node to weight
				decrease_weight(root,out,tr_list->weight);
				//keys[tr_list->dest]=tr_list->weight;
				// set the parent of this vertex as min vertex
				parent[tr_list->dest]=min->vertex;
				//visited[tr_list->dest]=true;
			}   
			tr_list = tr_list->next;
		}
		// avoid duplicate first entry
		if(total!=0){
			output.push_back(min->vertex);
			output.push_back(parent[min->vertex]);
		}
		// sum the cost of the tree so far with the min weight edge obtained.
		cost+=min->weight;
		total++;
	}
	// end timer.
	Time = clock() - Start; // time in micro seconds
	cout<<cost<<endl;
	while(!output.empty()){
		int s1=output.back();
		output.pop_back();
		int s2=output.back();
		output.pop_back();
		cout<<s1<<" "<<s2<<endl;
	}
	return Time;
}
// Prims algorithm implementation using a simple array to update weights.
// At each iteration we search the array to find the vertex associated with
// minimum weight and add it to the spanning tree.
clock_t prims_simple(struct Graph* g){
	clock_t Start, Time;
	// start the timer.
	Start = clock();
	int nodes=g->V;
	// array to store weights
	int keys[nodes];
	// array to mark parent of each vertex when visited
	int parent[nodes];
	// array to mark a vertex as visited so that we dont relook at the vertices.
	bool visited[nodes];
	list<int> output;
	// intialization of all the arrays
	for(int i=0;i<nodes;i++){
		keys[i]=2000;
		visited[i]=false;
		parent[i]=-1;
	}
	int cost=0,count=0;
	// assign weight of vertex 0 as 0 and start the algorithm.
	keys[0]=0;
	// while there are unexplored nodes
	while(count<nodes){
		int min=2000;
		int mark_node=0;
		// scan the keys array to get the vertex with minimum weight
		for(int i=0;i<nodes;i++){
			// if its weight is less than current weight and its not explored then
			// update min and remember vertex
			if(keys[i]<min && visited[i]!=true){
				min=keys[i];
				mark_node=i;
			}
		}
		// mark the minimum vertex as visited
		visited[mark_node]=true;
		// traverse the adjacency list of this vertex
		struct AdjListNode* tr_list = g->array[mark_node].head;
		while (tr_list){
			// if the vertex in the adj list is not explored and its weight is more than edge
			// weight then update it and set its parent.
			if(visited[tr_list->dest]==false and keys[tr_list->dest]>tr_list->weight){
				keys[tr_list->dest]=tr_list->weight;
				parent[tr_list->dest]=mark_node;
				//visited[tr_list->dest]=true;
			}   
			tr_list = tr_list->next;
		}
		// avoid duplicate first entry
		//pair<int, int> current_edge(mark_node,parent[mark_node]);
		if(count>0){
		// add the vertex and its parent in the spanning tree.
		//cout<<mark_node<<"-------"<<parent[mark_node]<<endl;
			output.push_back(mark_node);
			output.push_back(parent[mark_node]);
		}
		// sum the cost with new minimum edge we got.
		cost+=min;
		count++;
	}
	// end timer.
	Time = clock() - Start; // time in micro seconds
	//cout<<"end time : "<<clock()<<endl;
	//output.push_back(cost);
	cout<<cost<<endl;
	//output.pop_back();
	while(!output.empty()){
		int s1=output.back();
		output.pop_back();
		int s2=output.back();
		output.pop_back();
		cout<<s1<<" "<<s2<<endl;
	}
	return Time;
	//cout<<"cost "<<cost<<endl;
}
int main(int argc, char* argv[]) {
	if(argc==1){
		cout<<"Please specify the mode and paramenters to run the program!"<<endl;
		return 0;
	}
	if(strcmp(argv[1],"-s")==0 || strcmp(argv[1],"-f")==0){
		if(argv[2]!=NULL){
			string line;
			// read the input file
			ifstream fs (argv[2]);
			if (fs.is_open()){
				// get line by line data
				getline(fs,line);
				if(line!=""){
					// if file has some valid data, split the line and get integer values
					// here we get vertices number and edges number
					int num_vertices=atoi(line.substr(0,line.find(" ")).c_str());
					line.erase(0,line.find(" ")+1);
					int edges=atoi(line.c_str());
					//cout<<"indices and edges are "<<num_vertices<<" "<<edges<<endl;
					int edgeinput[edges*3];
					int line_count=0;
					// scan through file and populate edgeinput array with i,j,w values sequentially
					while ( getline (fs,line) ){
						string s = line;
						string delimiter = " ";
						size_t pos = 0;
						string token;
						while ((pos = s.find(delimiter)) != string::npos) {
							token = s.substr(0, pos);
							edgeinput[line_count++]=atoi(token.c_str());
							s.erase(0, pos + delimiter.length());
						}
						edgeinput[line_count++]=atoi(s.c_str());
					}
					// create graph with given vertices
					struct Graph* graph = createGraph(num_vertices);
					for(int i=2;i<edges*3;i+=3){
						// scan through the array and add the triplet as parameters to addEdge
						addEdge(graph,edgeinput[i-2],edgeinput[i-1],edgeinput[i]);
						//cout<<edgeinput[i-2]<<" "<<edgeinput[i-1]<<" "<<edgeinput[i]<<endl;
					}
					// if option is -s call simple scheme
					if(strcmp(argv[1],"-s")==0){
						//check for connectivity through dfs
						if(checkConnected(graph))
							clock_t time=prims_simple(graph);
							//cout<<"time for simple prims : "<<prims_simple(graph)<<endl;
						else
							cout<<"input graph is not connected!"<<endl;
					}
					// if option is -f call fibonacci scheme
					else if(strcmp(argv[1],"-f")==0){
						//check connectivity through dfs
						if(checkConnected(graph))
							clock_t time=prims_simple(graph);
							//cout<<"time for fibonacci prims : "<<prims_fscheme(graph)<<endl;
						else
							cout<<"input graph is not connected!"<<endl;
					}
				}
				//printGraph(graph);
				fs.close();
			}
			else{ 
				cout << "Unable to open file"<<endl;
			}
		}
	}
	else if(strcmp(argv[1],"-r")==0){
		int nodes,density;
		nodes=atoi(argv[2]);
		density=atoi(argv[3]);
		struct Graph* graph=generateRandomGraph(nodes,density);
		//printGraph(graph);
		//ans=checkConnected(graph);
		while(!checkConnected(graph)){
			//cout<<"not connected yet!";
			graph=generateRandomGraph(nodes,density);
		}
		//call prims simple scheme
		cout<<"time for simple prims : "<<prims_simple(graph)<<endl;
		cout<<"time for fibonacci prims : "<<prims_fscheme(graph)<<endl;
		//call fibinacci scheme
	}
	else{
		cout<<"Please specify the mode and paramenters to run the program!"<<endl;
	}
	//for(int i = 0; i < argc; i++)
	//	cout << "argv[" << i << "] = " << argv[i] << endl;
	return 0;
} 

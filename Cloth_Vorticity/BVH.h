#ifndef BVH_H_
#define BVH_H_
#include "BBox.h"
#include "BVHPrimitve.h"
#include "TIMER.h"
#include <Types.h>
#include <vector>
#include <queue>
#include <stdint.h>
#include <memory>
#include <iostream>

//! Node descriptor for the flattened tree
struct BVHFlatNode {
	BBox bbox;	//bounding box of this node
	int rightOffset; //the offset of its right child with regard to itself; the offset of the left child is always 1; If it is zero, then it is a leaf node

	// underlying Primitives if bbox overlaps
	int start;	//Primitive index in the Primitive list
	int nPrims; //(end - start): number of Primitives in the list covered by this node
};

struct BVHBuildEntry {
	// If not 0xfffffffc then this is the index of the parent in the flattened tree
	int parent;
	// The range of Primitives in the Primitive list covered by this node.
	int start, end;
};



template<typename Primitive>
class BVH {
public:

	BVH(std::vector<Primitive *> * Primitives, int leafSize = 1)
		: nNodes(0), nLeafs(0), leafSize(leafSize), build_prims(Primitives), flatTree(nullptr) {
		TIMER t;
		build();
		double constructionTime = t.Get_Time();
		//std::cout << "building time: " << constructionTime << "s"<< std::endl;
		buildBottomUpTree();
	}

	~BVH(){
		//delete[] flatTree;
	}
	
	void update() {
		TIMER t;
		for (size_t i = 0; i < bottomUpTree.size(); ++i) {
			int idx = bottomUpTree[i];
			BVHFlatNode& node = flatTree[idx];
			if (node.rightOffset == 0) {
				node.bbox = (*build_prims)[node.start]->getBBox();
				for (uint32_t i = 1; i < node.nPrims; ++i) {
					node.bbox.expandToInclude((*build_prims)[node.start + i]->getBBox());
				}
			}
			else {
				node.bbox = flatTree[idx + 1].bbox;
				node.bbox.expandToInclude(flatTree[idx + node.rightOffset].bbox);
			}
		}
		double updateTime = t.Get_Time();
		//std::cout << "updateTime : " << updateTime << "s" << std::endl;
	}

	bool testIntersection(const BBox& testBox, std::vector<Primitive> &candidates) const{
		const int stackSize = 1024;
		int stackptr = 0;
		int32_t todo[stackSize];
		candidates.clear();

		// "Push" on the root node to the working set
		todo[stackptr] = 0;
		while (stackptr >= 0) {
			// Pop off the next node to work on.
			int ni = todo[stackptr];
			stackptr--;
			const BVHFlatNode &node((flatTree)[ni]);

			// Is leaf -> Intersect
			if (node.rightOffset == 0) {
				for (int i = 0; i < node.nPrims; ++i) {
					Primitive* obj = (*build_prims)[node.start + i];
					if (testBox.overlap(obj->getBBox())){
						candidates.push_back(*obj);
					}
				}
			}
			else{ // Not a leaf
				bool hitc0 = testBox.overlap((flatTree)[ni + 1].bbox);
				bool hitc1 = testBox.overlap((flatTree)[ni + node.rightOffset].bbox);
				if (hitc0)
					todo[++stackptr] = ni + 1;
				if (hitc1)
					todo[++stackptr] = ni + node.rightOffset;
			}
			assert(stackptr < stackSize);
		}

		if (!candidates.empty())
			return true;
		else
			return false;
	}
	Primitive* getPrimitive(int i) const{
		assert(i< (*build_prims).size());
		return (*build_prims)[i];
	}
	BVHFlatNode getFlatTreeNode(int i) const{
		assert(i<nNodes);
		return flatTree[i];
	}
	int getTreeSize(){
		return nNodes;
	}

private:
	//! Build the BVH tree out of build_prims, it is a depth-first traversal
	void build(){
		const int stackSize = 1024;
		int stackptr = 0;
		BVHBuildEntry todo[stackSize];
		const int Untouched = 0xffffffff;
		const int TouchedTwice = 0xfffffffd;

		// Push the root
		todo[stackptr].start = 0;
		todo[stackptr].end = build_prims->size();
		todo[stackptr].parent = 0xfffffffc;
		stackptr++;

		BVHFlatNode node;
		std::vector<BVHFlatNode> buildnodes;
		buildnodes.reserve(build_prims->size() * 2);
		parentTable.reserve(build_prims->size() * 2);

		while (stackptr > 0) {
			// Pop the next item off of the stack
			BVHBuildEntry &bnode(todo[--stackptr]);
			int start = bnode.start;
			int end = bnode.end;
			int nPrims = end - start;

			nNodes++;
			node.start = start;
			node.nPrims = nPrims;
			node.rightOffset = Untouched;

			// Calculate the bounding box for this node
			BBox bb((*build_prims)[start]->getBBox());
			BBox bc((*build_prims)[start]->getCentroid());
			for (int p = start + 1; p < end; ++p) {
				bb.expandToInclude((*build_prims)[p]->getBBox());
				bc.expandToInclude((*build_prims)[p]->getCentroid());
			}
			node.bbox = bb;

			// If the number of primitives at this point is less than the leaf
			// size, then this will become a leaf. (Signified by rightOffset == 0)
			if (nPrims <= leafSize) {
				node.rightOffset = 0;
				nLeafs++;
			}

			buildnodes.push_back(node);

			// Child touches parent...
			// Special case: Don't do this for the root.
			if (bnode.parent != 0xfffffffc) {
				buildnodes[bnode.parent].rightOffset--;

				// When this is the second touch, this is the right child.
				// The right child sets up the offset for the flat tree.
				if (buildnodes[bnode.parent].rightOffset == TouchedTwice) {
					buildnodes[bnode.parent].rightOffset = nNodes - 1 - bnode.parent;
				}
				parentTable.push_back(bnode.parent);
			}else{
				parentTable.push_back(-1);
			}

			// If this is a leaf, no need to subdivide.
			if (node.rightOffset == 0)
				continue;

			// Set the split dimensions
			int split_dim = bc.maxDimension();

			// Split on the center of the longest axis
			float split_coord = .5f * (bc.min[split_dim] + bc.max[split_dim]);

			// Partition the list of Primitives on this split
			int mid = start;
			for (int i = start; i<end; ++i) {
				if ((*build_prims)[i]->getCentroid()[split_dim] < split_coord) {
					std::swap((*build_prims)[i], (*build_prims)[mid]);
					++mid;
				}
			}

			// If we get a bad split, just choose the center...
			if (mid == start || mid == end) {
				mid = start + (end - start) / 2;
			}

			// Push right child
			todo[stackptr].start = mid;
			todo[stackptr].end = end;
			todo[stackptr].parent = nNodes - 1;
			stackptr++;

			// Push left child
			todo[stackptr].start = start;
			todo[stackptr].end = mid;
			todo[stackptr].parent = nNodes - 1;
			stackptr++;

			assert(stackptr < stackSize);
		}

		flatTree = new BVHFlatNode[nNodes];
		// Copy the temp node data to a flat array
		for (int n = 0; n < nNodes; ++n)
			flatTree[n] = buildnodes[n];
	}

	void buildBottomUpTree() {
		bottomUpTree.clear();
		bottomUpTree.reserve(nNodes);

		std::queue<int> queque;
		std::vector<int> hitTimes(nNodes, 0);
		for (int n = 0; n < nNodes; ++n) {
			BVHFlatNode& flatNode = flatTree[n];
			if (flatNode.rightOffset == 0) {
				bottomUpTree.push_back(n);
				uint32_t parent = parentTable[n];

				++hitTimes[parent];
				if (hitTimes[parent] == 2) {
					queque.push(parent);
				}
			}
		}

		while (!queque.empty()) {
			int n = queque.front();
			queque.pop();

			bottomUpTree.push_back(n);
			int parent = parentTable[n];
			if (parent != -1) {
				++hitTimes[parent];
				if (hitTimes[parent] == 2) {
					queque.push(parent);
				}
			}
		}
	}

private:
	// build primitives
	std::vector<Primitive *> *build_prims;
	// Fast Traversal System
	BVHFlatNode *flatTree;

	int nNodes, nLeafs, leafSize;

	//for update
	std::vector<int> bottomUpTree;
	std::vector<int> parentTable;
};







//shared_ptr version

//! A Bounding Volume Hierarchy system 
//! requirement: Primitive must implement two methods: getBBox(), getCentroid()

/*
template<typename Primitive>
class BVH {
public:

	BVH(){}

	BVH(std::shared_ptr<std::vector<std::shared_ptr<Primitive> > > Primitives, int leafSize = 1)
		: nNodes(0), nLeafs(0), leafSize(leafSize), build_prims(Primitives), flatTree(NULL) {
		build();
	}

	~BVH(){}

	bool testIntersection(const BBox& testBox, std::vector<Primitive> &candidates) const{
		const int stackSize = 1024;
		int stackptr = 0;
		int32_t todo[stackSize];
		candidates.clear();

		// "Push" on the root node to the working set
		todo[stackptr] = 0;
		while (stackptr >= 0) {
			// Pop off the next node to work on.
			int ni = todo[stackptr];
			stackptr--;
			const BVHFlatNode &node((*flatTree)[ni]);

			// Is leaf -> Intersect
			if (node.rightOffset == 0) {
				for (int i = 0; i < node.nPrims; ++i) {
					std::shared_ptr<Primitive> obj = (*build_prims)[node.start + i];
					if (testBox.overlap(obj->getBBox())){
						candidates.push_back(*obj);
					}
				}
			}else{ // Not a leaf
				bool hitc0 = testBox.overlap((*flatTree)[ni + 1].bbox);
				bool hitc1 = testBox.overlap((*flatTree)[ni + node.rightOffset].bbox);
				if (hitc0) 
					todo[++stackptr] = ni + 1;
				if (hitc1) 
					todo[++stackptr] = ni + node.rightOffset;
			}
			assert(stackptr < stackSize);
		}

		if (!candidates.empty())
			return true;
		else
			return false;
	}
	std::shared_ptr<Primitive> getPrimitive(int i) const{
		return (*build_prims)[i];
	}
	BVHFlatNode getFlatTreeNode(int i) const{
		return (*flatTree)[i];
	}
	int getTreeSize(){
		return nNodes;
	}

private:
	//! Build the BVH tree out of build_prims, it is a depth-first traversal
	void build(){
		const int stackSize = 1024;
		int stackptr = 0;
		BVHBuildEntry todo[stackSize];
		const int Untouched = 0xffffffff;
		const int TouchedTwice = 0xfffffffd;

		// Push the root
		todo[stackptr].start = 0;
		todo[stackptr].end = build_prims->size();
		todo[stackptr].parent = 0xfffffffc;
		stackptr++;

		BVHFlatNode node;
		std::vector<BVHFlatNode> buildnodes;
		buildnodes.reserve(build_prims->size() * 2);

		while (stackptr > 0) {
			// Pop the next item off of the stack
			BVHBuildEntry &bnode(todo[--stackptr]);
			int start = bnode.start;
			int end = bnode.end;
			int nPrims = end - start;

			nNodes++;
			node.start = start;
			node.nPrims = nPrims;
			node.rightOffset = Untouched;

			// Calculate the bounding box for this node
			BBox bb( (*build_prims)[start]->getBBox());
			BBox bc( (*build_prims)[start]->getCentroid());
			for (int p = start + 1; p < end; ++p) {
				bb.expandToInclude((*build_prims)[p]->getBBox());
				bc.expandToInclude((*build_prims)[p]->getCentroid());
			}
			node.bbox = bb;

			// If the number of primitives at this point is less than the leaf
			// size, then this will become a leaf. (Signified by rightOffset == 0)
			if (nPrims <= leafSize) {
				node.rightOffset = 0;
				nLeafs++;
			}

			buildnodes.push_back(node);

			// Child touches parent...
			// Special case: Don't do this for the root.
			if (bnode.parent != 0xfffffffc) {
				buildnodes[bnode.parent].rightOffset--;

				// When this is the second touch, this is the right child.
				// The right child sets up the offset for the flat tree.
				if (buildnodes[bnode.parent].rightOffset == TouchedTwice) {
					buildnodes[bnode.parent].rightOffset = nNodes - 1 - bnode.parent;
				}
			}

			// If this is a leaf, no need to subdivide.
			if (node.rightOffset == 0)
				continue;

			// Set the split dimensions
			int split_dim = bc.maxDimension();

			// Split on the center of the longest axis
			float split_coord = .5f * (bc.min[split_dim] + bc.max[split_dim]);

			// Partition the list of Primitives on this split
			int mid = start;
			for (int i = start; i<end; ++i) {
				if ((*build_prims)[i]->getCentroid()[split_dim] < split_coord) {
					std::swap((*build_prims)[i], (*build_prims)[mid]);
					++mid;
				}
			}

			// If we get a bad split, just choose the center...
			if (mid == start || mid == end) {
				mid = start + (end - start) / 2;
			}

			// Push right child
			todo[stackptr].start = mid;
			todo[stackptr].end = end;
			todo[stackptr].parent = nNodes - 1;
			stackptr++;

			// Push left child
			todo[stackptr].start = start;
			todo[stackptr].end = mid;
			todo[stackptr].parent = nNodes - 1;
			stackptr++;

			assert(stackptr < stackSize);
		}

		flatTree = std::make_shared<std::vector<BVHFlatNode> >(std::vector<BVHFlatNode>(nNodes));
		// Copy the temp node data to a flat array
		for (int n = 0; n < nNodes; ++n)
			(*flatTree)[n] = buildnodes[n];
	}

private:
	// build primitives
	std::shared_ptr<std::vector<std::shared_ptr<Primitive> > > build_prims;
	// Fast Traversal System
	std::shared_ptr<std::vector<BVHFlatNode> > flatTree;

	int nNodes, nLeafs, leafSize;
};


//create BVH from mesh
static BVH<BVHPointPrimitive> createBVHFromPointMesh(float *points, int numPoints){
	std::vector<std::shared_ptr<BVHPointPrimitive> > primitives;
	for (size_t i = 0; i < numPoints; i++){
		primitives.push_back(std::make_shared<BVHPointPrimitive>(BVHPointPrimitive(points, i)));
	}
	BVH<BVHPointPrimitive> bvh(std::make_shared<std::vector<std::shared_ptr<BVHPointPrimitive> > >(primitives));
	return bvh;
}

static BVH<BVHTrianglePrimitive> creatBVHFromTriangleMesh(float *points, int *triIndex, int numTriangles){
	std::vector<std::shared_ptr<BVHTrianglePrimitive> > primitives;
	for (size_t i = 0; i < numTriangles; i++){
		primitives.push_back(std::make_shared<BVHTrianglePrimitive>(BVHTrianglePrimitive(points, triIndex[3 * i + 0], triIndex[3 * i + 1], triIndex[3 * i + 2])));
	}
	BVH<BVHTrianglePrimitive> bvh(std::make_shared<std::vector<std::shared_ptr<BVHTrianglePrimitive> > >(primitives));
	return bvh;
}

*/
#endif

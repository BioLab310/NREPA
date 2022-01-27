#ifndef _Bplustree_h_
#define _Bplustree_h_
#include "basic.h"

#define M 151
#define MID 76

struct Node
{
	int32_t	nodesize;
	uint32_t	data[M];
	char b[M];
	uint8_t		isleaf;
	struct Node *	child[M];
	struct Node *	brother;
	struct Node *	parent;
};
struct para_cmp
{
	uint32_t max_len;
	char * ref;
	uint32_t ref_len;
};
struct Node * NodeInitial();
void treeDestroy(struct Node * r);
uint32_t cmp(uint32_t a,\
		uint32_t b,\
		struct para_cmp para);
int32_t binarySearchPos(struct Node * p,\
		uint32_t val,\
		struct para_cmp para);
void leaf_shift_one(struct Node *p, int32_t i);
void node_split(struct Node ** root, \
		struct Node *p_father,\
		struct Node * p_current,\
		uint32_t val,\
		char b_insert,\
		struct para_cmp para);
void insert(struct Node ** root,\
		struct Node * tmp_p,\
		int val,\
		char b_insert,\
		struct para_cmp para);
struct Node* Find_Minimal_Node(struct Node* p_root);
#endif

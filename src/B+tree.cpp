#include "B+tree.h"

struct Node * NodeInitial()
{
	struct Node* p = (struct Node *)malloc(sizeof(struct Node));
	p->isleaf = 1;
	p->nodesize = 0;
	p->brother = NULL;
	p->parent = NULL;
	for (int i = 0; i<M; i++)
	{
		p->child[i] = NULL;
	}
	return p;
}
void treeDestroy(struct Node * r)
{
	if(r!=NULL)
	{
		for(int32_t i=0;i<r->nodesize;i++)
		{
			if(r->child[i]!=NULL)
			{
				treeDestroy(r->child[i]);
			}
		}
		free(r);
	}
}
uint32_t cmp(uint32_t a, uint32_t b,struct para_cmp para)
{
	char * p=para.ref;
	uint32_t len=para.ref_len;

	uint32_t x = 0;
	uint32_t i = 0;
	while (1)
	{
		if (p[a + i]<p[b + i])
		{
			x = 1;
			break;
		}
		if (p[a + i]>p[b + i])
		{
			break;
		}
		if ((a + i == (len - 1)) && (b + i == (len -1)))
		{
			return 2;
			break;
		}
		if (a + i == len - 1)
		{
			x = 1;
			break;
		}
		if (b + i == len - 1)
		{
			break;
		}
		if(para.max_len!=0)
		{
			if(i==para.max_len)
			{
				return 0;
			}
		}
		i++;
	}
	return x;
}
int32_t binarySearchPos(struct Node * p, uint32_t val, struct para_cmp para)
{
	if(p->nodesize==0)
	{		
		return 0;
	}
	else if(cmp(p->data[p->nodesize-1],val,para)==1)
	{
		return p->nodesize;
	}
	else if(cmp(val, p->data[0],para)==1)
	{
		return 0;
	}
	else
	{
		int32_t a=0;
		int32_t b=p->nodesize-1;
		int32_t mid;
		mid = (a+b)/2;
		while(1)
		{
			if(cmp(val,p->data[mid],para))
			{
				b=mid;
			}
			else
			{				
				a=mid;
			}
			mid=(a+b)/2;
			if(mid==a)
			{
				return mid+1;
			}
		}
	}
}
void leaf_shift_one(struct Node *p, int32_t i)
{
	for (int32_t j = p->nodesize; j > i; j--)
	{
		p->data[j] = p->data[j - 1];
		p->b[j]=p->b[j-1];
	}
}
void node_split(struct Node ** root, \
		struct Node *p_father,\
		struct Node * p_current,\
		uint32_t val,\
		char b_insert,\
		struct para_cmp para)
{
	if (p_father->nodesize != M)
	{
		cout << "Error: from calling node split!" << endl;
	}

	//1)
	int32_t insert_pos;
	insert_pos=binarySearchPos(p_father,val,para);

	//3)
	struct Node * r = NodeInitial();
	r->isleaf = p_father->isleaf;
	r->parent = p_father->parent;
	r->brother = p_father->brother;
	p_father->brother = r;

	if (insert_pos >= MID)
	{
		p_father->nodesize = MID;
		for (int32_t i = MID; i<insert_pos; i++)
		{
			r->data[i - MID] = p_father->data[i];
			r->b[i - MID] = p_father->b[i];
			r->child[i - MID] = p_father->child[i];
		}
		r->data[insert_pos - MID] = val;
		r->b[insert_pos - MID]=b_insert;
		r->child[insert_pos - MID] = p_current;
		for (int32_t i = insert_pos; i<M; i++)
		{
			r->data[i - MID + 1] = p_father->data[i];
			r->b[i - MID + 1] = p_father->b[i];
			r->child[i - MID + 1] = p_father->child[i];
		}
		r->nodesize = M + 1 - MID;
	}
	else
	{
		for (int32_t i = MID - 1; i<M; i++)
		{
			r->data[i - MID + 1] = p_father->data[i];
			r->b[i - MID + 1] = p_father->b[i];
			r->child[i - MID + 1] = p_father->child[i];
		}
		r->nodesize = M + 1 - MID;

		for (int32_t i = MID - 2; i >= insert_pos; i--)
		{
			p_father->data[i + 1] = p_father->data[i];
			p_father->b[i + 1] = p_father->b[i];
			p_father->child[i + 1] = p_father->child[i];
		}
		p_father->data[insert_pos] = val;
		p_father->b[insert_pos] = b_insert;
		p_father->child[insert_pos] = p_current;
		p_father->nodesize = MID;
	}

	//3.5
	for (int32_t i = 0; i<r->nodesize; i++)
	{
		if (r->child[i] != NULL)
		{
			r->child[i]->parent = r;
		}
	}

	//4)
	struct Node * new_root;
	if (p_father->parent == NULL)
	{
		new_root = NodeInitial();
		new_root->isleaf = 0;

		new_root->data[0] = p_father->data[p_father->nodesize - 1];
		new_root->child[0] = p_father;
		new_root->data[1] = r->data[r->nodesize - 1];
		new_root->child[1] = r;

		r->parent = new_root;
		p_father->parent = new_root;

		new_root->nodesize = 2;
		*root = new_root;
	}
	else
	{
		//2)
		int32_t parent_pos;
		parent_pos=binarySearchPos(p_father->parent,val,para);

		//5)
		p_father->parent->child[parent_pos] = r;
		uint32_t new_Mid_Bound;
		new_Mid_Bound = p_father->data[MID - 1];

		if (p_father->parent->nodesize == M)
		{
			node_split(root,p_father->parent, p_father, new_Mid_Bound,'0',para);
		}
		else
		{
			for (int32_t i = p_father->parent->nodesize-1; i>=parent_pos; i--)
			{
				p_father->parent->data[i + 1] = p_father->parent->data[i];
				p_father->parent->child[i + 1] = p_father->parent->child[i];
			}
			p_father->parent->data[parent_pos] = new_Mid_Bound;
			p_father->parent->child[parent_pos] = p_father;
			p_father->parent->nodesize++;
		}
	}
}
void insert(struct Node ** root,struct Node * tmp_p, int val,char b_insert,struct para_cmp para)
{
	if (tmp_p->isleaf)
	{
		if (tmp_p->nodesize == M)
		{
			node_split(root,tmp_p, NULL, val,b_insert,para);
		}
		else
		{
			int32_t pos;
			pos=binarySearchPos(tmp_p,val,para);
			if(pos==tmp_p->nodesize)
			{
				tmp_p->data[tmp_p->nodesize] = val;
				tmp_p->b[tmp_p->nodesize] = b_insert;
			}
			else
			{
				leaf_shift_one(tmp_p, pos);
				tmp_p->data[pos] = val;
				tmp_p->b[pos] = b_insert;
			}
			tmp_p->nodesize++;
		}
	}
	else
	{
		int pos;
		pos=binarySearchPos(tmp_p,val,para);

		if(pos==tmp_p->nodesize)
		{
			tmp_p->data[tmp_p->nodesize - 1] = val;
			insert(root,tmp_p->child[tmp_p->nodesize - 1], val,b_insert,para);
		}
		else
		{
			insert(root,tmp_p->child[pos], val,b_insert,para);
		}
	}
}
struct Node* Find_Minimal_Node(struct Node* p_root)
{
	struct Node* p_tmp=p_root;
	while(p_tmp->isleaf!=1)
	{
		p_tmp=p_tmp->child[0];
	}
	return p_tmp;
}

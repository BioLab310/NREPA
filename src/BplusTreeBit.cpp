/*
 * BplusTreeBit.cpp
 *
 *  Created on: Nov 12, 2018
 *      Author: bio
 */
#include "BplusTreeBit.h"

//#include "Binary_Search.h"

void cpyData(struct nodeBit* a,\
		struct nodeBit b,\
		struct bit256KmerPara para)
{
	if(a->hashValue==NULL)
	{
		a->hashValue=(uint64_t*)malloc(sizeof(uint64_t)*para.kmer64Len);
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			a->hashValue[j]=b.hashValue[j];
		}
	}
	else
	{
		for(uint32_t j=0;j<para.kmer64Len;j++)
		{
			a->hashValue[j]=b.hashValue[j];
		}
	}
	a->arrayID=b.arrayID;
}
void Node_initial_bit(struct NodeBit * p)
{
	p->Node_Size=0;
	p->leaf_label=1;
	p->numAsChild=-1;

	p->parent=NULL;
	p->brother=NULL;
	for(uint32_t i=0;i<M;i++)
	{
		p->p_child[i]=NULL;
		p->data[i].hashValue=NULL;
	}
}
void Insert_Value_bit(struct NodeBit** p_root,\
		struct nodeBit insert_data,\
		struct bit256KmerPara para)
{
	struct NodeBit*p_root_tmp= * p_root;
		//insert lead node
	while(1)
	{
		if(p_root_tmp->leaf_label==1)
		{
			break;
		}
		uint32_t child_label=0;
		for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				uint32_t gi;
				if(i==0)
				{
					gi=0;
				}
				else
				{
					gi=i-1;
				}
				p_root_tmp=p_root_tmp->p_child[gi];
				child_label=1;
				break;
			}
		}
		if(child_label==0)
		{
			p_root_tmp=p_root_tmp->p_child[p_root_tmp->Node_Size-1];
		}
	}


	uint32_t pos;
	if(p_root_tmp->Node_Size==0)
	{
		pos=0;
	}
	else
	{
		int32_t pos_label=0;
		for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==2)
			{
				pos=i;
				pos_label=-1;
				cout << "p_root_tmp->data[i].arrayID : " << p_root_tmp->data[i].arrayID << " " << "insert_data.arrayID : " << insert_data.arrayID << endl;
				getchar();
				break;
			}
			else if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				pos=i;
				pos_label=1;
				break;
			}
		}
		if(pos_label==-1)
		{
			cout << "Insert B+tree error!" << endl;
			return;
		}
		else if(pos_label==0)
		{
			pos=p_root_tmp->Node_Size;
		}
	}

	if(p_root_tmp->Node_Size==M)
	{
		Divide_Node_bit(p_root,p_root_tmp,insert_data,pos,NULL,0,0,para);
	}
	else
	{

		for(uint32_t i=pos;i<p_root_tmp->Node_Size;i++)
		{
			uint32_t l=p_root_tmp->Node_Size-1-(i-pos);
			cpyData(&(p_root_tmp->data[l+1]),p_root_tmp->data[l],para);
		}
		cpyData(&(p_root_tmp->data[pos]),insert_data,para);
		p_root_tmp->Node_Size++;

		if(pos==0&&p_root_tmp->parent!=NULL)
		{
			NodeBit* p_parent_update=p_root_tmp->parent;
			cpyData(&(p_parent_update->data[p_root_tmp->numAsChild]),p_root_tmp->data[0],para);
			NodeBit* p1=p_root_tmp;
			NodeBit* p2=p_parent_update;
			while(p1->numAsChild==0&&p2->parent!=NULL)
			{
				cpyData(&(p2->parent->data[p2->numAsChild]),p2->data[0],para);
				p1=p2;
				p2=p2->parent;
			}
		}
	}
}
void Insert_Node_bit(struct NodeBit** p_root,\
		struct NodeBit* p_inserted_node,\
		struct nodeBit insert_data,\
		struct NodeBit* p_middle,\
		struct bit256KmerPara para)
{
	struct NodeBit*p_root_tmp= p_inserted_node;

	uint32_t pos;
	uint32_t pos_label=0;
	for(uint32_t i=0;i<p_root_tmp->Node_Size;i++)
	{
		if(cmp256BitKmer(insert_data.hashValue,p_root_tmp->data[i].hashValue,para.kmer64Len)==0)
		{
			pos=i;
			pos_label=1;
			break;
		}
	}
	if(pos_label==0)
	{
		pos=p_root_tmp->Node_Size;
	}
	if(p_root_tmp->Node_Size==M)
	{
		Divide_Node_bit(p_root,p_root_tmp,insert_data,pos,p_middle,0,0,para);
	}
	else
	{
		for(uint32_t i=pos;i<p_root_tmp->Node_Size;i++)
		{
			uint32_t l=p_root_tmp->Node_Size-1-(i-pos);
			cpyData(&(p_root_tmp->data[l+1]),p_root_tmp->data[l],para);

			p_root_tmp->p_child[l+1]=p_root_tmp->p_child[l];
			p_root_tmp->p_child[l+1]->numAsChild++;
		}

		cpyData(&(p_root_tmp->data[pos]),insert_data,para);

		p_root_tmp->p_child[pos]=p_middle;
		p_middle->numAsChild=pos;
		p_middle->parent=p_root_tmp;

		p_root_tmp->Node_Size++;

		if(pos==0&&p_root_tmp->parent!=NULL)
		{
			NodeBit* p_parent_update=p_root_tmp->parent;
			cpyData(&(p_parent_update->data[p_root_tmp->numAsChild]),p_root_tmp->data[0],para);
			NodeBit* p1=p_root_tmp;
			NodeBit* p2=p_parent_update;
			while(p1->numAsChild==0&&p2->parent!=NULL)
			{
				cpyData(&(p2->parent->data[p2->numAsChild]),p2->data[0],para);
				p1=p2;
				p2=p2->parent;
			}
		}
	}
}
void Divide_Node_bit(struct NodeBit** p_root,\
		struct NodeBit* p_divide_node,\
		struct nodeBit insert_data,\
		uint32_t pos,\
		struct NodeBit* insert_child,\
		uint32_t ad_label,\
		uint32_t poslist_label,\
		struct bit256KmerPara para)
{
	struct NodeBit * parent_tmp;
	parent_tmp=p_divide_node->parent;

	struct NodeBit * p_add=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_add);
	p_add->leaf_label=p_divide_node->leaf_label;
	p_add->Node_Size=M/2+1;
	if(p_divide_node->leaf_label==1)
	{
		p_add->brother=p_divide_node->brother;
		p_divide_node->brother=p_add;
	}

	if(pos<=p_divide_node->Node_Size-M/2-1)
	{
		uint32_t shift=(p_divide_node->Node_Size-M/2-1);
		for(uint32_t i=shift;i<p_divide_node->Node_Size;i++)
		{
			cpyData(&(p_add->data[i-shift]),p_divide_node->data[i],para);
			p_add->p_child[i-shift]=p_divide_node->p_child[i];

			if(p_divide_node->p_child[i]!=NULL)
			{
				p_divide_node->p_child[i]->parent=p_add;
				p_divide_node->p_child[i]->numAsChild=i-shift;
			}
			p_divide_node->p_child[i]=NULL;
		}

		for(uint32_t i=pos;i<=shift-1;i++)
		{
			uint32_t l=shift-1-(i-pos);
			cpyData(&(p_divide_node->data[l+1]),p_divide_node->data[l],para);
			p_divide_node->p_child[l+1]=p_divide_node->p_child[l];

			if(p_divide_node->p_child[l+1]!=NULL)
			{
				p_divide_node->p_child[l+1]->numAsChild=l+1;
			}
		}

		cpyData(&(p_divide_node->data[pos]),insert_data,para);
		p_divide_node->p_child[pos]=insert_child;

		if(insert_child!=NULL)
		{
			insert_child->parent=p_divide_node;
			insert_child->numAsChild=pos;
		}
		p_divide_node->Node_Size=p_divide_node->Node_Size-M/2;

		if(pos==0&&p_divide_node->parent!=NULL)
		{
			NodeBit* p_parent_update=p_divide_node->parent;
			cpyData(&(p_parent_update->data[p_divide_node->numAsChild]),p_divide_node->data[0],para);
			NodeBit* p1=p_divide_node;
			NodeBit* p2=p_divide_node->parent;
			while(p1->numAsChild==0&&p2->parent!=NULL)
			{
				cpyData(&(p2->parent->data[p2->numAsChild]),p2->data[0],para);

				p1=p2;
				p2=p2->parent;
			}
		}

		if(parent_tmp==NULL)
		{
			update_root_bit(p_root,p_add,para);
		}
		else
		{
			Insert_Node_bit(p_root,parent_tmp,p_add->data[0],p_add,para);
		}
	}
	else
	{
		uint32_t shift=(p_divide_node->Node_Size-M/2);
		for(uint32_t i=p_divide_node->Node_Size-M/2;i<=pos-1;i++)
		{
			cpyData(&(p_add->data[i-shift]),p_divide_node->data[i],para);
			p_add->p_child[i-shift]=p_divide_node->p_child[i];

			if(p_divide_node->p_child[i]!=NULL)
			{
				p_divide_node->p_child[i]->parent=p_add;
				p_divide_node->p_child[i]->numAsChild=i-shift;
			}
			p_divide_node->p_child[i]=NULL;
		}

		cpyData(&(p_add->data[pos-shift]),insert_data,para);

		p_add->p_child[pos-shift]=insert_child;
		if(insert_child!=NULL)
		{
			insert_child->parent=p_add;
			insert_child->numAsChild=pos-shift;
		}
		for(uint32_t i=pos;i<=p_divide_node->Node_Size-1;i++)
		{
			cpyData(&(p_add->data[i-shift+1]),p_divide_node->data[i],para);

			p_add->p_child[i-shift+1]=p_divide_node->p_child[i];

			if(p_divide_node->p_child[i]!=NULL)
			{
				p_divide_node->p_child[i]->parent=p_add;
				p_divide_node->p_child[i]->numAsChild=i-shift+1;
			}
			p_divide_node->p_child[i]=NULL;
		}
		p_divide_node->Node_Size=p_divide_node->Node_Size-M/2;

		if(parent_tmp==NULL)
		{
			update_root_bit(p_root,p_add,para);
		}
		else
		{
			Insert_Node_bit(p_root,parent_tmp,p_add->data[0],p_add,para);
		}
	}
}
void update_root_bit(struct NodeBit** p_root,\
		struct NodeBit* p2,\
		struct bit256KmerPara para)
{
	struct NodeBit *root_tmp=*p_root;
	struct NodeBit * p_tmp=(struct NodeBit*)malloc(sizeof(struct NodeBit));
	Node_initial_bit(p_tmp);

	p_tmp->Node_Size=2;
	p_tmp->leaf_label=0;

	cpyData(&(p_tmp->data[0]),root_tmp->data[0],para);
	p_tmp->p_child[0]=root_tmp;
	root_tmp->numAsChild=0;
	root_tmp->parent=p_tmp;

	cpyData(&(p_tmp->data[1]),p2->data[0],para);
	p_tmp->p_child[1]=p2;
	p2->numAsChild=1;
	p2->parent=p_tmp;

	*p_root=p_tmp;
}
void destory_tree_bit(struct NodeBit* p_root,\
		struct bit256KmerPara para)
{
	if(p_root!=NULL)
	{
		for(uint32_t i=0;i<p_root->Node_Size;i++)
		{
			destory_tree_bit(p_root->p_child[i],para);
			if(p_root->data[i].hashValue != NULL)
			{
				free(p_root->data[i].hashValue);
				p_root->data[i].hashValue = NULL;
			}
		}
		free(p_root);
	}
}
struct NodeBit* Find_Minimal_Node_bit(struct NodeBit* p_root)
{
	struct NodeBit* p_tmp=p_root;
	while(p_tmp->leaf_label!=1)
	{
		p_tmp=p_tmp->p_child[0];
	}
	return p_tmp;
}
int64_t MappingHashValueToID_bit(struct NodeBit* p_root,\
		uint64_t* hashValue,\
		struct bit256KmerPara para)
{
	struct NodeBit* p_tmp = p_root;
	while(p_tmp!=NULL&&p_tmp->leaf_label!=1)
	{
		uint64_t label=0;
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(i!=0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				p_tmp=p_tmp->p_child[i-1];
				label=1;
				break;
			}
			else if(i==0&&cmp256BitKmer(hashValue,p_tmp->data[i].hashValue,para.kmer64Len)==0)
			{
				return -1;
			}

		}
		if(label==0)
		{
			p_tmp=p_tmp->p_child[p_tmp->Node_Size-1];
		}
	}

	if(p_tmp==NULL)
	{
		return -1;
	}
	else
	{
		for(uint64_t i=0;i<p_tmp->Node_Size;i++)
		{
			if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==1)
			{
				return -1;
			}
			else if(cmp256BitKmer(p_tmp->data[i].hashValue,hashValue,para.kmer64Len)==2)
			{
				return p_tmp->data[i].arrayID;
			}
		}
		return -1;
	}
}


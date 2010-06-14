#ifndef _AVL_TREE_H_
#define _AVL_TREE_H_


#include <stdio.h>
//
// AVL Tree Implementation
// (c) Copyright 2004 David Gu.  All rights reserved.
//
//
// The template class AVL::AVLTreeNode (used in AVL::Tree) requires the following:
//
//   T::T ( T &);
//   T::~T ();
//   bool T::operator < ( T & rhs) ;
//   bool T::operator == ( T & rhs) ;
// 
// For printing purposes one must declare:
//
//   ostream & operator << ( ostream &,  T &);
//


#define	TREE_MAX_HEIGHT 64

//
// The basic unit of currency in a tree are the nodes that comprise it.
//
template <class T>
class __declspec( dllexport ) AVLTreeNode
{
public :

// This is were we keep the data we want to store in each node.
// It is  because if you change it while it is in the tree structure
// you compromise the integrity of the tree.
// It is public because the Tree class must have access to it in order
// to return it after being found with the found_node function.

     T * data;



// Each node has two children: left and right.  If they are both NULL then
// the node is a leaf node.  Otherwise, it's an interior node.

    AVLTreeNode<T> * left, * right;

private :

// The height is computed to be: 0 if NULL, 1 for leaf nodes, and the maximum
// height of the two children plus 1 for interior nodes.
// This is used to keep the tree balanced.

    int height;

    void compute_height ()
    {
        height = 0;
        if (left != NULL && left -> height > height)
            height = left -> height;
        if (right != NULL && right -> height > height)
            height = right -> height;
        height += 1;
    }

// The ructor is private because the nodes are self allocating.

    AVLTreeNode (T *  inData)
        : data (inData), left (NULL), right (NULL), height (1)
    {
    }

public :

// Recursively delete the children if this node is being nuked.

    ~AVLTreeNode ()
    {
        delete left;
        delete right;
    }

// Recursively insert some data into the tree then balance it on the way up.

    AVLTreeNode<T> * insert_node (T  *  inData)
    {
        if (this == NULL)
            return new AVLTreeNode<T> (inData);

        if ( *inData < *data)
            left = left -> insert_node (inData);
        else
            right = right -> insert_node (inData);
        return balance ();
    }

// Recursively find some data in the tree and if found return a pointer
// to the node containing the data.  If not found then return NULL.

    AVLTreeNode<T> * find_node (T * inData) 
    {
        if (this == NULL)
            return NULL;

        if ( *inData == *data)
            return this;

        if (*inData < *data)
            return left -> find_node (inData);
        else
            return right -> find_node (inData);
    }

// Recursively search the tree for some data and if found remove (delete) it.
// When you remove an interior node the right child must be place right of
// the right most child in the left sub-tree.
// Remember to balance the tree on the way up after removing a node.

    AVLTreeNode<T> * remove_node (T * inData)
    {
        if (this == NULL)
            return NULL;

        // we found the data we were looking for

        if ( *inData == *data)
        {
            // save the children

            AVLTreeNode<T> * tmp = left -> move_down_righthand_side (right);

            // by setting the children to NULL, we delete exactly one node.

            left = NULL;
            right = NULL;
            delete this;

            // return the reorganized children

            return tmp;
        }

        if ( *inData < *data)
            left = left -> remove_node (inData);
        else
            right = right -> remove_node (inData);
        return balance ();
    }


private :


// move_down_righthand_side is the remove_node helper function:
//
// Recursively find the right most child in a sub-tree and put
// the "rhs" sub-tree there.
// Re-balance the tree on the way up.

    AVLTreeNode<T> * move_down_righthand_side (AVLTreeNode<T> * rhs)
    {
        if (this == NULL)
            return rhs;

        right = right -> move_down_righthand_side (rhs);
        return balance ();
    }

//
// Balancing a tree (or sub-tree) requires the AVL algorithm.
//
// If the tree is out of balance left-left, we rotate the node to the right.
// If the tree is out of balance left-right, we rotate the left child to the
// left and then rotate the current node right.
// If the tree is out of balance right-left, we rotate the right child to the
// right and then rotate the current node left.
// if the tree is out of balance right-right, we rotate the node to the left.
//

    AVLTreeNode<T> * balance ()
    {
        int d = difference_in_height ();

        // only rotate if out of balance
        if (d < -1 || d > 1)
        {
            // too heavy on the right
            if (d < 0)
            {
                // if right child is too heavy on the left,
                // rotate right child to the right
                if (right -> difference_in_height () > 0)
                    right = right -> rotate_right ();

                // rotate current node to the left
                return rotate_left ();
            }
            // too heavy on the left
            else
            {
                // if left child is too heavy on the right,
                // rotate left child to the left
                if (left -> difference_in_height () < 0)
                    left = left -> rotate_left ();

                // rotate current node to the right
                return rotate_right ();
            }
        }

        // recompute the height of each node on the way up
        compute_height ();

        // otherwise, the node is balanced and we simply return it
        return this;
    }

// ** balancing helper functions **

    AVLTreeNode<T> * exchange_left (AVLTreeNode<T> * & r, AVLTreeNode<T> * node)
    {
        r = left;
        left = node -> balance ();
        return balance ();
    }

    AVLTreeNode<T> * exchange_right (AVLTreeNode<T> * & l, AVLTreeNode<T> * node)
    {
        l = right;
        right = node -> balance ();
        return balance ();
    }

    int difference_in_height ()
    {
        int left_height = (left != NULL) ? left -> height : 0;
        int right_height = (right != NULL) ? right -> height : 0;
        return left_height - right_height;
    }

    AVLTreeNode<T> * rotate_left ()
    {
        return right -> exchange_left (right, this);
    }

    AVLTreeNode<T> * rotate_right ()
    {
        return left -> exchange_right (left, this);
    }

};

//
// Cover class for maintaining the tree.
//
// Since AVLTreeNode<T> is self allocating and self deleting, the AVLTree<T> class
// ensures that only qualified calls are made.
//
// AVLTree<T> is the public interface to the AVL Tree code.
// AVLTreeNode<T> is not meant to be used by the public.
//
// This code makes use of the somewhat dubious practice of calling a member
// function with a NULL "this" pointer.  We will not run into problems since
// we have no virtual member functions in AVLTreeNode<T>.
//

template <class T>
class __declspec( dllexport ) AVLTree
{
private :

    AVLTreeNode<T> * root;

public :

    AVLTree ()
    {
        root = NULL;
    }

    ~AVLTree ()
    {
        delete root;
    }

    void insert (T * inData)
    {
        root = root -> insert_node (inData);
    }

    T * find (T * inData) 
    {
        AVLTreeNode<T> * found = root -> find_node (inData);
        if (found != NULL)
            return  found -> data;
        else
            return NULL;
    }

    void remove (T * inData)
    {
        root = root -> remove_node (inData);
    }

	AVLTreeNode<T> * Root() { return root; };	
};


// Used for traversing an AVL tree. 
template <class T> 
class __declspec( dllexport ) AVLTreeIterator
{
	public:
		AVLTreeIterator( AVLTree<T> & tree );
		~AVLTreeIterator(){};
		void operator++();
		T *  operator*();
		bool end() { return m_finished; };
		void reset();

	private:
		bool m_finished;                    //finished ?
		bool m_initialized;					//initialized ?
		int  m_top;							//top of stack
		AVLTreeNode<T> * m_pointer;				//current AVLTreeNode
		AVLTreeNode<T> * m_stack[TREE_MAX_HEIGHT];	// Descended trees. 
		AVLTreeNode<T> * m_root;                   // root of the tree
};

template<class T>
AVLTreeIterator<T>::AVLTreeIterator( AVLTree<T> & tree )
{
	m_initialized = true;
	m_finished    = false;

	m_top  = 0;
	m_root = tree.Root();
	m_pointer = m_root;
	

	for (;;)
    {
		/* T2. */
		while (m_pointer != NULL)
		{
			/* T3. */
			m_stack[m_top ++ ] = m_pointer;
			m_pointer = m_pointer->left;
		}
		
		/* T4. */
		if (m_top == 0)
		{
			m_finished = true;
			m_initialized = false;
			return;
		}

		m_pointer = m_stack[--m_top];
		return;
    }

};

template<class T>
void AVLTreeIterator<T>::reset()
{
	m_initialized = true;
	m_finished    = false;

	m_top = 0;
	m_pointer = m_root;
	

	for (;;)
    {
		/* T2. */
		while (m_pointer != NULL)
		{
			/* T3. */
			m_stack[m_top ++ ] = m_pointer;
			m_pointer = m_pointer->left;
		}
		
		/* T4. */
		if (m_top == 0)
		{
			m_finished = true;
			m_initialized = false;
			return;
		}

		m_pointer = m_stack[--m_top];
		return;
    }

};


template <class T> 
T * AVLTreeIterator<T>::operator *()
{
	return m_pointer->data;
}

template<class T>
void AVLTreeIterator<T>::operator++()
{	
	/* Uses Knuth's algorithm 2.3.1T (inorder Traverserersal). */
	/* T5. */
	m_pointer = m_pointer->right;
	
	for (;;)
    {
		/* T2. */
		while ( m_pointer != NULL)
		{
			/* T3. */
			m_stack[m_top ++ ] = m_pointer;
			m_pointer = m_pointer->left;
		}
		
		/* T4. */
		if (m_top == 0)
		{
			m_initialized = false;
			m_finished = true;
			return;
		}

		m_pointer = m_stack[--m_top];
		
		/* T5. */
		return;
    }
};



#endif

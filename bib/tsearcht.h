// tsearcht.h 3.1 97/08/05 01:57:49
//
//	AVL-Suchbaum


# if !defined( TSEARCHT_H_INCLUDED )

#define TSEARCHT_H_INCLUDED

#include <iostream.h>

template <class T>
class TreeNode
{
 public:

   TreeNode( T init );

   TreeNode<T> *&left( ) { return leftTree; }
   TreeNode<T> *&right( ) { return rightTree; }

   T &data( ) { return nodeData; }
   void setData( const T &to ) { nodeData = to; }
   unsigned count( ) { return nodeCount; }
   void incCount( ) { ++nodeCount; }
   void decCount( ) { --nodeCount; }

   int compare(const T &to);
   void set( const TreeNode<T> &to );

   void deleteTree( );
   unsigned height( ) const;
   void print( unsigned h, ostream &out );
   void foreach(void (* f)(void *, void *), void *data = 0);

   void rotateRR( );
   void rotateLL( );
   void rotateLRLR( );
   void rotateRLRL( );

 private:

   TreeNode<T> *leftTree;
   TreeNode<T> *rightTree;

   unsigned nodeCount;
   T nodeData;
};

template <class T>
class AVLNode : public TreeNode<T>
{
 public:

   AVLNode( T init );

   AVLNode<T> *&left( )  { return (AVLNode<T> *&) TreeNode<T>::left( ); }
   AVLNode<T> *&right( ) { return (AVLNode<T> *&) TreeNode<T>::right( ); }

   int balance( ) { return nodeBalance; }
   void setBalance( int to = 0 ) { nodeBalance = to; }

 protected:

  int nodeBalance;
};

template <class T>
class SearchTree
{
   typedef TreeNode<T> Node;

 public:

   SearchTree( );
   virtual ~SearchTree( );

   unsigned height( ) const;
   

   virtual void insert( const T &elem ) { ins( elem, tree ); }
   virtual void remove( const T &elem ) { rm( elem, tree ); }
   
   T * find(const T &like) const;
   void foreach(void (* f)(void *, void *), void *data = 0);

   friend ostream &operator<<( ostream &out, const SearchTree<T> &t )
      { if ( t.tree ) t.tree->print( 0, out ); return out; }

 protected:

   TreeNode<T> *tree;

 private:

   static void ins( const T &x, TreeNode<T> *&node );
   static void del( TreeNode<T> *&node, TreeNode<T> *&branch );
   static void rm( const T &x, TreeNode<T> *&node );
};

template <class T>
class AVLTree : public SearchTree<T>
{
  typedef AVLNode<T> Node;

 public:

   AVLTree() : SearchTree<T>() {}

   virtual void insert( const T &elem );
   virtual void remove( const T &elem );

 protected:

   static int isAVL( AVLNode<T> *&node );

   static void rebalanceLeft( AVLNode<T> *&node );
   static void rebalanceRight( AVLNode<T> *&node );
   static void rebalance( AVLNode<T> *&node, int &change, int sign );
   static void ins( const T &x, AVLNode<T> *&node, int &change );

   static void balanceLeft( AVLNode<T> *&node, int &change );
   static void balanceRight( AVLNode<T> *&node, int &change );
   static void del( AVLNode<T> *&node, AVLNode<T> *&branch, int &change );
   static void rmNode( AVLNode<T> *&node, int &change );
   static void rm( const T &x, AVLNode<T> *&node, int &change );
};

#if defined(TEMPL_INCL_DEF)
#   if defined(CC_SOURCE)
#      include <tsearcht.cc>
#   else
#      include <tsearcht.cpp>
#   endif
#endif

#endif









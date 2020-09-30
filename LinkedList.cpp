#include "stdafx.h"

node * LinkedList::first() const
{
	return _first;
}

LinkedList::LinkedList()
{
	_first = NULL;
}

LinkedList::LinkedList(float x)
{
	_first = new node(x);
}

LinkedList::LinkedList(node * newnode)
{
	_first = newnode;
}

void LinkedList::destroy() {
	node *current = _first;
	while (current) {
		node *temp = current;
		current = current->_next;
		delete temp;
	}
	_first = NULL;
}

LinkedList::~LinkedList(void)
{
	destroy();
}

void LinkedList::insertAtFront(float x) {
	node* n = new node(x);
	n->_next = _first;
	_first = n;
}

void LinkedList::insertAtFront(node *n) {
	n->_next = _first;
	_first = n;
}

void LinkedList::Print()const {
	node* current = _first;
	while (current) {
		cout << current->_data;
		if (current->_next)
			cout << " -> ";
		current = current->_next;
	}
	cout << endl;
}

int LinkedList::size()const {
	int l = 0;
	node* current = _first;
	while (current) {
		l++;
		current = current->_next;
	}
	return l;
}

node* LinkedList::find(float item) const {
	node* current = _first;
	while (current) {
		if (current->_data == item)
			return current;
		current = current->_next;
	}
	return NULL;
}

void LinkedList::insertAtEnd(float x){
	node* n = new node(x);
	if (_first == NULL) {
		_first = n;
		return;
	}
	node* current = _first;
	while (current->_next) 
		current = current->_next;
	current->_next = n;
}

void LinkedList::insertAtEnd(node * n){
	n->_next = NULL;
	if (_first == NULL) {
		_first = n;
		return;
	}
	node* current = _first;
	while (current->_next)
		current = current->_next;
	current->_next = n;
}

void LinkedList::insertAfter(float x,  node * des){
	node* n = new node(x);
	n->_next = des->_next;
	des->_next = n;
}

void LinkedList::insertAfter(node * n,  node * des){
	n->_next = des->_next;
	des->_next = n;
}

void LinkedList::deleteFromFront()
{
	if (_first == NULL)
		return;
	node* temp = _first;
	_first = _first->_next;
	delete temp;
}

void LinkedList::deleteFromEnd(){
	if (_first == NULL)  //List is Empty
		return;
	if (_first->_next == NULL) { //List has One Node
		delete _first;
		_first = NULL;
		return;
	}
	node* secondLast = _first;
	node* currentPrev = NULL;
	while (secondLast->_next->_next) {
		secondLast = secondLast->_next;
	}
	delete secondLast->_next;
	secondLast->_next = NULL;
}

void LinkedList::deleteNode(node * n) {
	if (_first == NULL || n == NULL)  //List is Empty
		return;
	if (_first == n) {
		_first = _first->_next;
		delete n;
		return;
	}
	node* current = _first;
	while (current->_next != n) {
		current = current->_next;
	}
	current->_next = n->_next;
	delete n;
}

void LinkedList::deleteNode(float item)
{
	node* n = find(item);
	deleteNode(n);
}


void LinkedList::sort(){
	float t;
	for (node* icurrent = _first; icurrent->_next != NULL; icurrent = icurrent->_next)
		for (node* jcurrent = icurrent->_next; jcurrent != NULL; jcurrent = jcurrent->_next)
			if (icurrent->_data < jcurrent->_data)
				swap1(icurrent->_data , jcurrent->_data);
				//swap2(icurrent,jcurrent);

}


void LinkedList::swap2(node* xnode,node* ynode){
   node *prevX = NULL, *currX = _first;
   while (currX!= xnode){
	   prevX=currX;
	   currX=currX->_next;
   }
   node* prevY = NULL, *currY = _first;
   while (currY!= ynode){
	   prevY = currY;
	   currY=currY->_next;
   }
   
   if (currX == NULL ||currY==NULL)
	   return;
   
   if (prevX != NULL)
	   prevX->_next = currY;
   else
	   _first = currY;
   
   if (prevY != NULL)
	   prevY->_next = currX;
   else
	   _first = currX;
   
   node * temp= currY->_next;
   currY->_next = currX->_next;
   currX->_next=temp;
}




void LinkedList::swap1(float& x,float& y){
	float t=x;
	x=y;
	y=t;
}
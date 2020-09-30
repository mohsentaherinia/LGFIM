#include<iostream>
#include<conio.h>
using namespace std;
#pragma once
class LinkedList;
class node{
	friend class LinkedList;
private:
	float _data;
	node* _next;
public:
	node::node(float x) {
		_data = x;
		_next = NULL;
	}
	node::node(float x, node* link) {
		_data = x;
		_next = link;
	}
	float getData()const{
		return _data;
	}
	node* getNext()const {
		return _next;
	}
	void setData(const float &data) {
		_data = data;
	}
	void setNext(const float &next) {
		_data = next;
	}
};


class LinkedList
{
private:
	node* _first;
	void swap1(float&,float&);
	void swap2(node* icurrent,node* jcurrent);
public:
	node* first()const;
	LinkedList();
	LinkedList(float x);
	LinkedList(node* node);
	void destroy();
	~LinkedList(void);
	void insertAtFront(float x);
	void insertAtFront(node * n);
	void Print()const;
	int size()const;
	void insertAtEnd(float x);
	void insertAtEnd(node * n);
	void insertAfter(float x,  node *des);
	void insertAfter(node* x,  node *des);
	void deleteFromFront();
	void deleteFromEnd();
	void deleteNode(node * n);
	void deleteNode(float item);
	node* find(float item)const;
	void sort();
};


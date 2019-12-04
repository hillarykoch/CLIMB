#ifndef _LLIST_H
#define _LLIST_H

#include "mynode.h"

using namespace lemon;

class llist {
    private:
        node *head;
        node *tail;

    public:
        llist() {
            head = NULL;
            tail = NULL;
        }

        void insert_first(ListDigraph::Node value) {
            node *temp = new node;
            temp->data = value;
        
            if(head == NULL) {
                head = temp;
                tail = temp;
                head->next = NULL;
                tail->next = NULL;
                temp = NULL;
            } else {
                temp->next = head;  // point the new node to head
                head = temp; // move head back one, to the new node
            }
        }
        
        void insert_last(ListDigraph::Node value) {//int value) {
            node *temp = new node;
            temp->data = value;
            temp->next = NULL; // temp doesnt point to anything
            if(head == NULL) {
                head = temp;
                tail = temp;
                temp = NULL;
            } else {
                tail->next = temp; // point current tail to temp
                tail = temp; // move tail to the new node
            }        
        }
        
        void insert_at(int index, ListDigraph::Node value) {
            // currently no exception handling for trying to access an invalid index
            node *temp = new node;
            node *curr = new node;
            node *prev = new node;
            temp->data = value;
        
            curr = head->next;
            prev = head;
        
            int count = 0;
            while(count < index - 1) {
                prev = curr;
                curr = curr->next;
                count++;   
            }
            prev->next = temp;
            temp->next = curr;
        }
        
        void remove_first() {
            if(head == tail) {
                // If there is only one node, point head and tail to null
                head = NULL;
                tail = NULL;
            }
            
            // currently not handling exception when list is empty
        
            node *temp = new node;
            temp = head;
            head = head->next;
            delete temp;
        }
        
        void remove_last() {
            if (head == tail) {
                // If there is only one node, point head and tail to null
                std::cout << "shouldn't be here" << std::endl;
                head = NULL;
                tail = NULL;
            }
        
            // currently not handling exception when list is empty
            node *curr = new node;
            //node *prev = new node;
            curr = head;
            while(curr->next != NULL) { // start at head and iterate over llist
              //prev = curr;
              tail=curr;
              curr = curr->next;    
            }
            //tail=prev;
            //prev->next=NULL;
            tail->next=NULL;
            delete curr;
        }
        
        void remove_at(int index) {
            // currently no exception handling for trying to access an invalid index
            node *curr = new node;
            node *foll = new node;
        
            curr = head;
            int count = 0;
            while(count < index - 1) {
                foll = foll->next;
                curr = curr->next;
                count++;   
            }
            curr->next = foll->next;
            delete foll;
        }
        
        const ListDigraph::Node& operator[](int index) const {
            // currently not handling exception when index is too big
            node *curr = new node;
            curr = head;
            int count = 0;
            while(count < index) {
                // start at head and iterate over llist
                // this indexing starts at 0!!!
                count++;
                curr = curr->next;    
            }
            return curr->data;
        }
        
        const int len() const {
            if(head == NULL) { // if there is nothing
                return 0;
            } else {
                node *curr = new node;
                curr = head;
                int count = 1;
                while(curr->next != NULL) {
                    count++;
                    curr = curr->next;
                }
                return count;
            }
        }
};


#endif
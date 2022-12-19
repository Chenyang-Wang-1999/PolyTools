/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-19
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#ifndef LINKLIST_TOOLS_HPP
#define LINKLIST_TOOLS_HPP

template<typename NodeType>
void destroy_tree(NodeType * node_ptr)
{
    if(node_ptr == NULL)
    {
        return;
    }
    else
    {
        if(node_ptr->next == NULL)
        {
            node_ptr->~NodeType();
        }
        else
        {
            destroy_tree(node_ptr->next);
            node_ptr->~NodeType();
        }
    }
}

template<typename NodeType>
NodeType * pop_next_node(NodeType * curr_node_ptr)
{
    NodeType * next_node_ptr = curr_node_ptr -> next;
    assert(next_node_ptr != NULL);
    curr_node_ptr -> next = next_node_ptr -> next;
    next_node_ptr -> next = NULL;
    return next_node_ptr;
}

template<typename NodeType>
void remove_next_node(NodeType * curr_node_ptr)
{
    NodeType * next_node_ptr = pop_next_node<NodeType>(curr_node_ptr);
}

template<typename NodeType>
void insert_after_ptr(NodeType * curr_node_ptr, NodeType * new_node)
{
    new_node -> next = curr_node_ptr -> next;
    curr_node_ptr->next = new_node;
}

template<typename NodeType>
void traverse_from_node(NodeType * head_ptr, void (*funcall)(NodeType *))
{
    NodeType * curr_ptr = head_ptr;
    while(curr_ptr != NULL)
    {
        funcall(curr_ptr);
        curr_ptr = curr_ptr->next;
    }
}

template<typename NodeType, typename DataType>
void traverse_from_node(NodeType * head_ptr, void (*funcall)(NodeType *, DataType *), DataType* fundata)
{
    NodeType * curr_ptr = head_ptr;
    while(curr_ptr != NULL)
    {
        funcall(curr_ptr, fundata);
        curr_ptr = curr_ptr -> next;
    }
}

#endif
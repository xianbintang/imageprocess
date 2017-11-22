#include "List.h"

/*需要头结点*/
void delete_node(List L, CvRect item);

static List find_min(List L);

List create_list()
{
	List L;
	L = (struct Node *)malloc(sizeof(Node));

	if (L == NULL) {
		fprintf(stderr, "Error! out of space!\n");
	}

	L->next = NULL;

	return L;
}

void push_back(List L, CvRect i)
{
	struct Node * tmp = (struct Node *)malloc(sizeof(Node));
	if (tmp == NULL) {
		fprintf(stderr, "Error! out of space!\n");
	}


	tmp->item = i;
	L->next = tmp;
	tmp->next = NULL;
}

/*连头结点也算在内了*/
int count_node(List L)
{
	int count = 0;
	while (L != NULL) {
		count++;
		L = L->next;
	}

	return count;
}

bool is_empty(List L)
{
	return (L->next == NULL);
}

/*清空一个链表*/
void empty_list(List L)
{
	L->next = NULL;
}

List sort(List L)
{
	List sorted_list = create_list();
	List head = sorted_list;
	List min_node = NULL;

	while (L->next != NULL) {
		min_node = find_min(L);
		push_back(sorted_list, min_node->item);
		delete_node(L, min_node->item);
		sorted_list = sorted_list->next;
	}

#if 0
	while (head->next != NULL) {
		printf("x: %d\n", head->next->item.x);
		head = head->next;
	}
#endif

	return head;
}

/*找到L下一个开始的最小节点*/
List find_min(List L)
{
	if (L == NULL) {
		fprintf(stderr, "L is NULL\n");
		exit(-1);
	}
	int min = 100000;
	List min_node = NULL;
	
	while (L->next != NULL) {
		if (L->next->item.x < min) {
			min = L->next->item.x;
			min_node = L->next;
		}
		L = L->next;
	}
	return min_node;
}

/*头结点一起传进来才行*/
void delete_node(List L, CvRect item)
{
	if (L == NULL)
		return;

	while (L->next != NULL) {
		if (L->next->item.x == item.x) {
			L->next = L->next->next;
				return;
		}

		L = L->next;
	}
}
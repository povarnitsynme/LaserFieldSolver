#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <cctype>
#include <cstdlib>
using namespace std;


enum node_t{UNKNOWNNODE, OPERATOR, NUMBER, ID, FUNC1, FUNC2};
enum oper_t{PLUS, MINUS, MULTIPL, DIVISION, NUM_OPERATORS};

const int precedence[NUM_OPERATORS] = {-5,-5,-4,-4};

typedef double (*fn1)(double); //function with 1 argument
typedef double (*fn2)(double, double); //function with 2 arguments


struct Node
{

	Node* m_parent;
	Node* m_left;
	Node* m_right;
	bool m_sign; //0 if no sign or plus, 1 if there is minus

	typedef double number_t;
	
	node_t m_type;
	void* m_content;

	Node()
	{
		m_type = UNKNOWNNODE; 
		m_content = (void*)0;
		m_parent = m_left = m_right = 0;
		m_sign = false;
	}
	
	Node(const oper_t op)
	{
		m_type = OPERATOR;
		m_content = new oper_t(op);
		m_parent = m_left = m_right = 0;
		m_sign = false;
	}
	Node(const number_t& number)
	{
		m_type = NUMBER;
		m_content = new number_t(number);
		m_parent = m_left = m_right = 0;
		m_sign = false;
	}
	Node(fn1 func)
	{
		m_type = FUNC1;
		m_content = reinterpret_cast<void*>(func);
		m_parent = m_left = m_right = 0;
		m_sign = false;
	}
	Node(fn2 func)
	{
		m_type = FUNC2;
		m_content = reinterpret_cast<void*>(func);
		m_parent = m_left = m_right = 0;
		m_sign = false;
	}

	Node(const string& id)
	{
		m_type = ID;
		m_content = new string(id);
		m_parent = m_left = m_right = 0;
		m_sign = false;
	}

	~Node()
	{
		switch(m_type)
		{
		case OPERATOR: delete (oper_t*)m_content; break;
		case NUMBER:   delete (number_t*)m_content; break;
		case ID:       delete (string*)m_content; break;
		case FUNC1:
		case FUNC2:	break; //no need to delete a function
		}
	}
};


class Expression
{
public:
	static map<string, fn1> funcs1;
	static map<string, fn2> funcs2;

private:
	char* m_expr_str;
	const char* m_pcur;
	Node* m_root;
	
	void eat_spaces()
	{
		if(m_pcur)
			while(isspace(*m_pcur)) ++m_pcur;
	}

	void deleteSubtree(const Node* node)
	{
		if(node)
		{
			deleteSubtree(node->m_left);
			deleteSubtree(node->m_right);
			delete node;
		}
	}
		
	double calculateSubtree(const Node* const node, const double& time, int* err_code = NULL) const
	{
		err_code?*err_code=0:0;

		if(!node) { err_code?*err_code=1:0; return 777.777;}
		switch(node->m_type)
		{
		case NUMBER: 
			return *(double*)(node->m_content) * ((node->m_sign) ? -1 : 1);;
		case ID: 
			return time * ((node->m_sign) ? -1 : 1);
		case OPERATOR:
			{
				int err_code1=0, err_code2=0;
				double a = calculateSubtree(node->m_left,  time, &err_code1);
				double b = calculateSubtree(node->m_right, time, &err_code2);
				if(err_code1 + err_code2)
				{
					err_code?*err_code=3:0;
					return -111.111;
				}
				switch(*(oper_t*)(node->m_content))
				{
				case PLUS:
					return (a+b) * ((node->m_sign) ? -1 : 1);
				case MINUS:
					return (a-b) * ((node->m_sign) ? -1 : 1);
				case MULTIPL:
					return (a*b) * ((node->m_sign) ? -1 : 1);
				case DIVISION:
					return (a/b) * ((node->m_sign) ? -1 : 1);
				}
			}break;
		case FUNC1:
			{
				double d = calculateSubtree(node->m_right, time, err_code);
				double res = ((fn1)(node->m_content))(d);
				return res * ((node->m_sign) ? -1 : 1);
			}
		case FUNC2:
			{
				double d1 = calculateSubtree(node->m_left, time, err_code);
				if(err_code && *err_code) return 0;
				double d2 = calculateSubtree(node->m_right, time, err_code);
				if(err_code && *err_code) return 0;
				double res = ((fn2)(node->m_content))(d1, d2);
				return res * ((node->m_sign) ? -1 : 1);
			}
		default:
			err_code?*err_code=4:0;
			return -888.888;
		}
		
		return -888.888;
	}

	Node* read_op()
	{
		eat_spaces();

		Node* new_node = (Node*)0;
		switch(*m_pcur)
		{
		case '+': ++m_pcur; new_node = new Node(PLUS);     break;
		case '-': ++m_pcur; new_node = new Node(MINUS);    break;
		case '*': ++m_pcur; new_node = new Node(MULTIPL);  break;
		case '/': ++m_pcur; new_node = new Node(DIVISION); break;
		}
		return new_node;
	}

	Node* read_id() //identifier contains [a-z0-9_], starts with a letter and doesn't contain '_'
	{
		eat_spaces();
		string id;
		if(isalpha(*m_pcur) || *m_pcur == '_')
		{
			while(isalpha(*m_pcur) || isdigit(*m_pcur) || *m_pcur == '_')
				id += *m_pcur++;

			return new Node(id);
		}
		return (Node*)0;
	}
	Node* read_number()
	{
		eat_spaces();
		Node* new_node = (Node*)0;
		char* tail;
		double d = strtod(m_pcur, &tail);
		if(tail != m_pcur)
		{
			m_pcur = tail;
			new_node = new Node(d);
		}
		return new_node;
	}

	bool read_c(char c)
	{
		eat_spaces();
		if(c == *m_pcur)
		{ 
			++m_pcur;
			return true;
		}
		else
			return false;
	}

	bool read_func1name(fn1& f)
	{
		eat_spaces();
		const char* p = m_pcur;
		string funcname;

		if( *p == '_' || isalpha(*p) )
			funcname += *p++;
		else
			return false;
		
		while( *p == '_' || isalnum(*p) )
			funcname += *p++;

		map<string, fn1>::const_iterator it = funcs1.find(funcname);
		if(it != funcs1.end())
		{
			m_pcur = p;
			f = it->second;
			return true;
		}
		return false;
	}

	bool read_func2name(fn2& f)
	{
		eat_spaces();
		const char* p = m_pcur;
		string funcname;

		if( *p == '_' || isalpha(*p) )
			funcname += *p++;
		else
			return false;
		
		while( *p == '_' || isalnum(*p) )
			funcname += *p++;

		map<string, fn2>::const_iterator it = funcs2.find(funcname);
		if(it != funcs2.end())
		{
			m_pcur = p;
			f = it->second;
			return true;
		}
		return false;

	}

	Node* read_func1()
	{
		eat_spaces();
		
		const char* p_save = m_pcur;

		Node* argnode  = (Node*)0;
		Node* funcnode = (Node*)0;
		fn1 f;

		if(read_func1name(f) && read_c('(') && (argnode = read_tree()) && read_c(')'))
		{
			funcnode = new Node(f);
			funcnode->m_left = (Node*)0;
			funcnode->m_right = argnode;
			argnode->m_parent = funcnode;
			return funcnode;
		}
		else
		{
			m_pcur = p_save;
			delete argnode;
			return (Node*)0;
		}
	}

	Node* read_func2()
	{
		eat_spaces();

		const char* p_save = m_pcur;

		Node* arg1node  = (Node*)0;
		Node* arg2node  = (Node*)0;
		Node* funcnode  = (Node*)0;
		fn2 f;

		if(read_func2name(f) && read_c('(') && (arg1node = read_tree()) &&
			read_c(',') && (arg2node = read_tree()) && read_c(')'))
		{
			funcnode = new Node(f);
			funcnode->m_left  = arg1node;
			funcnode->m_right = arg2node;
			arg1node->m_parent = funcnode;
			arg2node->m_parent = funcnode;
			return funcnode;
		}
		else
		{
			m_pcur = p_save;
			delete arg1node;
			delete arg2node;
			return (Node*)0;
		}
	}

	Node* read_node()
	{
		eat_spaces();
		
		Node* new_node = (Node*)0;

		if( new_node = read_number() )
			return new_node;
		else if( new_node = read_func1() )
			return new_node;
		else if( new_node = read_func2() )
			return new_node;
		else if( new_node = read_id() )
		{
			if(*(string*)(new_node->m_content) == "time") //only variable "time" admitted
			{
				return new_node;
			}
			else
			{
				delete new_node;
				return (Node*)0;
			}
		}
		else if( read_c('(') && (new_node = read_tree()) && read_c(')') )
			return new_node;
		else
		{
			delete new_node;
			return (Node*)0;
		}
	}

	Node* read_tree()
	{
		eat_spaces();
		Node* rmnode, *root;

		bool sign = false;
		if( read_c('-') ) sign = true;
		else if( read_c('+') ) {}

		rmnode = root = read_node();
		if(sign)
			root->m_sign = true;

		if(!root) return (Node*)0;

		Node* nextop, *nextnode;
		while(nextop = read_op())
		{
			nextnode = read_node();
			if(nextnode)
			{
				nextnode->m_parent = nextop;
				nextop->m_right = nextnode;

				Node* where_to_insert = rmnode->m_parent; //parent is always op, not number
				while( where_to_insert &&
							precedence[*(oper_t*)(nextop->m_content)] <= 
							precedence[*(oper_t*)(where_to_insert->m_content)] )	
					where_to_insert = where_to_insert->m_parent;
				
				if(where_to_insert == (Node*)0) //came to the top
				{
					nextop->m_left = root;
					root->m_parent = nextop;
					root = nextop;
				}
				else
				{
					where_to_insert->m_right->m_parent = nextop;
					nextop->m_left = where_to_insert->m_right;
					nextop->m_parent = where_to_insert;
					where_to_insert->m_right = nextop;
				}
				rmnode = nextnode;
				
			}
			else 
			{
				deleteSubtree(root);
				delete nextop;
				return (Node*)0;
			}
		}
		return root;
	}

private:
		static double min(const double a, const double b)
		{
			return (a<b)? a:b;
		}
		static double max(const double a, const double b)
		{
			return (a>b)? a:b;
		}
public:
	Expression(string expr)
	{
		size_t len = expr.length();
		m_expr_str = new char[len+1];
		memcpy(m_expr_str, expr.c_str(), len*sizeof(char));
		m_expr_str[len] = '\0';
		m_pcur = m_expr_str;
		m_root = (Node*)0;

		if(funcs1.empty())
		{
			funcs1["sin"] = &::sin;
			funcs1["cos"] = &::cos;
			funcs1["tan"] = &::tan;
			funcs1["atan"] = &::atan;
			funcs1["abs"] = &::fabs;
			funcs1["exp"] = &::exp;
			funcs1["log10"] = &::log10;
			funcs1["log"]  = &::log;
			funcs1["sqrt"]= &::sqrt;


			//funcs2.insert(make_pair("pow", mypow));
			funcs2.insert(make_pair("min", &Expression::min));
			funcs2.insert(make_pair("max", &Expression::max));

		}

		if( !build_tree() )
		{
			deleteSubtree(m_root);
			m_root = 0;
		}

	}

	~Expression() 
	{
		delete[] m_expr_str;
		deleteSubtree(m_root);
	}

	double calculate(const double& time, int* err_code = NULL) const
	{
		return calculateSubtree(m_root, time, err_code);
	}

	
	bool build_tree()
	{
		m_pcur = m_expr_str;
		m_root = read_tree();
		if(*m_pcur == 0) return true; // come to the end
		else return false;            // there may be odd symbols at the end
	}

	bool isCorrect() const
	{
		return m_root != (Node*)0;
	}

};

map<string, fn1> Expression::funcs1 = map<string, fn1>();
map<string, fn2> Expression::funcs2 = map<string, fn2>();

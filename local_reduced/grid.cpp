#include "grid.h"

Grid::Grid()
{
	root = NULL;
	size = 0;	
}

void Grid::add_node(vector<double> v)
{
	if(root == NULL)
	{
		root = new Node;
		root->set_points(v);
	}
	else
	{
		Node* parent = root;
		Node* child = root->get_next();
		while(child != NULL)
		{
			parent = parent->get_next();
			child = child->get_next();
		}
		child = new Node;
		child->set_points(v);
		parent->set_next( child );
	}
	size++;
}

vector< vector<double> > Grid::build_grid()
{
	if( size < 2 )
	{
		Node* temp;
		vector< vector<double> > full_grid;
		vector<double> vec,vec2;

		temp = root;
		vec = temp->get_points();

		for(int i = 0 ; i < vec.size() ; i++)
		{
			vec2.clear();
			vec2.push_back( vec[i] );
			full_grid.push_back(vec2);
		}

		return full_grid;
	}
	else
	{
		Node* temp;
		vector< vector<double> > full_grid;
		vector<double> vec,vec2;

		temp = root;
		vec = temp->get_points();

		for(int i = 0 ; i < vec.size() ; i++)
		{
			vec2.clear();
			vec2.push_back( vec[i] );
			full_grid.push_back(vec2);
		}

		while( temp->get_next() != NULL )
		{
			temp = temp->get_next();

			vec.clear();
			vec = temp->get_points();
			full_grid = kron_product(full_grid,vec);
		}
		return full_grid;
	}
}

vector< vector<double> > Grid::kron_product(vector< vector<double> > mat , vector<double> vec)
{
	vector< vector<double> > res_mat;
	vector<double> res_vec;
	vector<double> temp;
	
	for(int i = 0 ; i < mat.size() ; i++)
	{
		for(int j = 0 ; j < vec.size() ; j++)
		{
			temp.clear();
			temp = mat[i];
			for(int k = 0 ; k < temp.size() ; k++)
				temp[k] = temp[k] * 1;
			res_mat.push_back(temp);
		}
	}
	
	for(int i = 0 ; i < mat.size() ; i++)
	{
		for(int j = 0 ; j < vec.size() ; j++)
		{
			res_vec.push_back( vec[j] );
		}
	}

	for(int i = 0 ; i < res_mat.size() ; i++)
		res_mat[i].push_back( res_vec[i] );

	return res_mat;
}

int Grid::get_size()
{
	return size;
}

void Grid::print_grid()
{
	Node* temp = root;
	vector<double> points;
	int iter = 0;
	while(temp)
	{
		cout << "Node " << iter << " : " << endl;
		points = temp->get_points();
		for(int i = 0 ; i < points.size() ; i++)
			cout << points[i] << endl;
		temp = temp->get_next();
		iter++;
	}
}

void Grid::print_vec(vector<double> vec)
{
	cout << "----------------" << endl;
	for(int i = 0 ; i < vec.size() ; i++)
	{
		cout << vec[i] << " ";
	}
	cout << endl;
}

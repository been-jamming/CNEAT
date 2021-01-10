typedef struct gene gene;

struct gene{
	int from;
	int to;
	double weight;
	long int innovation_num;
	unsigned char activated;
};

typedef struct genome genome;

struct genome{
	gene *genes;
	int num_genes;
	int num_neurons;
	int species;
	double fitness;
};

typedef struct organism organism;

struct organism{
	double fitness;
	neuron **neurons;
	int num_neurons;
};

typedef struct species species;

struct species{
	int num_genomes;
	genome **genomes;
};

extern double add_connection_chance;
extern double mutate_chance;
extern double split_chance;
extern double toggle_activation_chance;
extern double random_weight_magnitude;

extern long int global_innovation_num;

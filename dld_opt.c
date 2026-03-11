#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

unsigned int xorshift_state;

void xorshift_seed(unsigned int seed) {
    if (seed == 0) {
        seed = time(NULL) ^ getpid();
    }
    xorshift_state = seed;
}

unsigned int xorshift32() {
    unsigned int x = xorshift_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    xorshift_state = x;
    return x;
}

double rand_double() {
    return xorshift32()/4294967296.0;
}

#define L 1048576 // Ancho del sustrato
#define Lm 3000 // Altura máxima de la caja de simulación
#define T 80 // Pasos de tiempo 
#define N 10 // Número de muestras


#define diameter 2  
#define SYSTEM_WIDTH (L * diameter) 
#define NUM_PARTICLES (L + T * L)

#define D_MAX 160          
#define L_MIN 1.0         
#define COLLISION_THRESHOLD 4

#define L_PARAM 4
#define CELL_SIZE L_PARAM
#define GRID_DIM_X (SYSTEM_WIDTH/CELL_SIZE)
#define GRID_DIM_Y (Lm/CELL_SIZE)
#define BLOCK_SIZE (L_PARAM*L_PARAM)
#define TOTAL_CELLS (GRID_DIM_X*GRID_DIM_Y)

#define MIN_LARGE_TREE_SIZE 100

#define LAST_PARTICLES_TO_SAVE L/20
double active_h_vector[LAST_PARTICLES_TO_SAVE + 1];

// Configuración de los tiempos a guardar

#define NUM_SNAPSHOTS 8
int snapshot_times[NUM_SNAPSHOTS] = {10, 20, 30, 40, 50, 60, 70, 80};

int job_id;

int W_grid[GRID_DIM_X][GRID_DIM_Y];
int next_k_value;
int Nk_counts[GRID_DIM_X*GRID_DIM_Y + 1];

int F_indices[(size_t)TOTAL_CELLS*BLOCK_SIZE];

int Y_grid[SYSTEM_WIDTH][Lm]; 
unsigned short int Omega_grid[SYSTEM_WIDTH][Lm];
unsigned short int Vicinity_grid[2*D_MAX + 1][2*D_MAX + 1];

// Variables globales

int particle_count, TIME[T];
long double M1[T], M2[T], M3[T], M4[T], M1_active[T], M2_active[T], M3_active[T], M4_active[T], results[30][2], 
total_results[NUM_SNAPSHOTS][30][2];
long double N_PARTICLES[T], total_rg_per_layer[T]; // Almacenar el número total de partículas promedio por capa
double particles_list[NUM_PARTICLES][2], interface_heights[L];
double collision_dist_sq;

// Array para acumular el RMS Thickness promedio por capa temporal

long double total_rms_per_layer[T];

// Variables globales para el análisis de árboles

#define MAX_TREES 400
#define MAX_N_BIN 1000000 // Máximo N_tree (masa) a promediar. Un árbol de masa mayor será ignorado

int umbral_minimo = 400;

// Arrays para promedio sobre N muestras

long double H_sum_binned[NUM_SNAPSHOTS][MAX_N_BIN];
int Count_binned[NUM_SNAPSHOTS][MAX_N_BIN];

#define MAX_GRID_CELLS 300000000 

int *grid_particles[MAX_GRID_CELLS];
int grid_counts[MAX_GRID_CELLS];      
int grid_capacities[MAX_GRID_CELLS];

// Spatial grid para búsqueda de vecinos

int grid_width;
int grid_height;
double grid_cell_size;

// BFS

int bfs_visited[NUM_PARTICLES];
int bfs_queue[NUM_PARTICLES];

void create_spatial_grid(double cell_size, double max_height) {

    grid_cell_size = cell_size;
    grid_width = (int)(SYSTEM_WIDTH/cell_size) + 1;
    grid_height = (int)(max_height/cell_size) + 1;
    
    memset(grid_particles, 0, sizeof(grid_particles));
    memset(grid_counts, 0, sizeof(grid_counts));
    memset(grid_capacities, 0, sizeof(grid_capacities));
}

void free_spatial_grid() {

    int total_cells = grid_width*grid_height;
    for (int i = 0; i < total_cells; i++) {
        if (grid_particles[i]) {
            free(grid_particles[i]);
            grid_particles[i] = NULL;
        }
    }
}

void add_to_spatial_grid(int particle_index) {

    double x = particles_list[particle_index][0];
    double y = particles_list[particle_index][1];
    
    int cell_x = (int)(x/grid_cell_size);
    int cell_y = (int)(y/grid_cell_size);
    
    if (cell_x < 0 || cell_x >= grid_width || 
        cell_y < 0 || cell_y >= grid_height) return;
    
    int cell_index = cell_y*grid_width + cell_x;

    if (cell_index >= MAX_GRID_CELLS) {
        fprintf(stderr, "Grid desbordada. cell_index=%d > MAX=%d\n", cell_index, MAX_GRID_CELLS);
        exit(1); 
    }
    
    if (grid_counts[cell_index] >= grid_capacities[cell_index]) {

        int new_capacity = (grid_capacities[cell_index] == 0) ? 8 : grid_capacities[cell_index]*2;

        grid_particles[cell_index] = realloc(grid_particles[cell_index], new_capacity*sizeof(int));
        grid_capacities[cell_index] = new_capacity;
    }
    
    grid_particles[cell_index][grid_counts[cell_index]++] = particle_index;
}

// Verificación de conexión

static inline int is_connected_fast(int p1_idx, int p2_idx, double diameter_sq) {

    double x1 = particles_list[p1_idx][0];
    double y1 = particles_list[p1_idx][1];
    double x2 = particles_list[p2_idx][0];
    double y2 = particles_list[p2_idx][1];
    
    double dy = y1 - y2;
    double dx = x1 - x2;
    
    if (dx > SYSTEM_WIDTH/2.0) dx -= SYSTEM_WIDTH;
    else if (dx < -SYSTEM_WIDTH/2.0) dx += SYSTEM_WIDTH;
    
    double dist_sq = dx*dx + dy*dy;
    
    return (dist_sq <= diameter_sq + 1e-9);
}

// Búsqueda de vecinos

void find_neighbors_fast(int particle_idx, double diameter_sq, int *neighbors, int *neighbor_count) {

    double x = particles_list[particle_idx][0];
    double y = particles_list[particle_idx][1];
    
    int cell_x = (int)(x/grid_cell_size);
    int cell_y = (int)(y/grid_cell_size);
    
    *neighbor_count = 0;
    
    for (int dy = -1; dy <= 1; dy++) {

        for (int dx = -1; dx <= 1; dx++) {

            int check_x = cell_x + dx;
            int check_y = cell_y + dy;
            
            if (check_x < 0) check_x += grid_width;
            if (check_x >= grid_width) check_x -= grid_width;
            
            if (check_y < 0 || check_y >= grid_height) continue;
            
            int cell_index = check_y * grid_width + check_x;
            
            for (int i = 0; i < grid_counts[cell_index]; i++) {

                int candidate = grid_particles[cell_index][i];
                
                if (candidate != particle_idx && 
                    is_connected_fast(particle_idx, candidate, diameter_sq)) {
                    neighbors[(*neighbor_count)++] = candidate;
                }
            }
        }
    }
}

// Análisis de los árboles

void analyze_trees_bfs(int total_particles, double diam, int sample_id, int current_time, int snap_idx) {

    double diameter_sq = diam*diam;

    // Abrir archivo de visualización

    char all_trees[100];

    sprintf(all_trees, "all_trees_raw_t=%d_run=%d.dat", current_time, job_id);

    FILE *f_raw = fopen(all_trees, "a");
    
    // Si el archivo está vacío (es la primera vez), escribimos la cabecera

    if (f_raw && ftell(f_raw) == 0) fprintf(f_raw, "# Sample_ID N H\n");

    memset(bfs_visited, 0, sizeof(int)*total_particles);

    // Encontrar altura máxima

    double max_height = 0;

    for (int i = 0; i < total_particles; i++) {

        if (particles_list[i][1] > max_height) max_height = particles_list[i][1];
    }
    
    // Spatial Grid 

    create_spatial_grid(diam*2.5, max_height + 10);

    for (int i = 0; i < total_particles; i++) add_to_spatial_grid(i);
    
    int neighbors[1000];
    int neighbor_count;

    // BFS Loop
    
    for (int i = L; i < total_particles; i++) {
        
        // Si ya visitamos esta partícula en un árbol previo, saltar

        if (bfs_visited[i]) continue;

        // Nuevo arbol

        int head = 0;
        int tail = 0;
        
        int current_tree_size = 0;
        double current_tree_max_h = 0;
        int touches_substrate = 0;

        // Iniciar cola con la partícula actual

        bfs_queue[tail++] = i;
        bfs_visited[i] = 1;

        while(head < tail) {

            int curr_idx = bfs_queue[head++]; // Dequeue
            
            // Procesar partícula

            current_tree_size++;
            double h = particles_list[curr_idx][1];
            if (h > current_tree_max_h) current_tree_max_h = h;

            // Buscar vecinos

            find_neighbors_fast(curr_idx, diameter_sq, neighbors, &neighbor_count);

            for (int j = 0; j < neighbor_count; j++) {

                int nbr_idx = neighbors[j];

                if (nbr_idx < L) touches_substrate = 1; // El sustrato no se añade a la cola, solo marca la condición

                else if (!bfs_visited[nbr_idx]) {

                    bfs_visited[nbr_idx] = 1; // Marcar
                    bfs_queue[tail++] = nbr_idx; // Enqueue

                }
            }
        }

        // Fin arbol
        // Guardar estadísticas si cumple condiciones

        if (touches_substrate && current_tree_size > umbral_minimo) {

            if (current_tree_size < MAX_N_BIN) {
                H_sum_binned[snap_idx][current_tree_size] += (long double)current_tree_max_h;
                Count_binned[snap_idx][current_tree_size]++;
            }

            // Datos crudos

            if (f_raw) fprintf(f_raw, "%d %.6f\n", current_tree_size, current_tree_max_h);
        }

    }

    if (f_raw) fclose(f_raw);
    free_spatial_grid();
}

// Guardar la media de las alturas para cada arbol de masa N

void save_averaged_tree_data() {
    
    // Bucle para guardar cada snapshot por separado

    for (int s = 0; s < NUM_SNAPSHOTS; s++) {
        
        int time_val = snapshot_times[s];
        char arch_tree[100];

        sprintf(arch_tree, "tree_scaling_avg_t=%d_L=%d_run=%d.dat", time_val, L, job_id);

        FILE *file_tree = fopen(arch_tree, "w");
        fprintf(file_tree, "# N_tree H_tree_avg Count\n");
        
        for (int n_tree = umbral_minimo; n_tree < MAX_N_BIN; n_tree++) {

            // Usamos el índice [s]

            if (Count_binned[s][n_tree] > 0) {

                long double avg_H = H_sum_binned[s][n_tree]/(long double)Count_binned[s][n_tree];
                fprintf(file_tree, "%d %.9Lf %d\n", n_tree, avg_H, Count_binned[s][n_tree]);
            }
        }

        fclose(file_tree);
    }
}

// Inicializar vectores de arboles

void init_tree_analysis() {

    // Limpiamos todo el bloque de memoria para los 4 tiempos

    memset(H_sum_binned, 0, sizeof(H_sum_binned));
    memset(Count_binned, 0, sizeof(Count_binned));
}

// Función para mostrar barra de progreso

void mostrar_barra_progreso(int actual, int total, const char* prefijo, time_t tiempo_inicio) {

    int porcentaje = (actual * 100)/total;

    // Calcular tiempo transcurrido
    time_t tiempo_actual = time(NULL);
    int segundos = (int)difftime(tiempo_actual, tiempo_inicio);
    int minutos = segundos/60;
    int segundos_restantes = segundos%60;

    // Limpiar la línea anterior y sobrescribir
    printf("\r\033[K");

    // Imprimir barra de progreso con tiempo en formato mm:ss
    printf("%s: [", prefijo);
    for (int j = 0; j < 60; j++) {
        if (j < porcentaje/2) printf("#");
        else printf(" ");
    }
    printf("] %3d%% | Tiempo: %02d:%02d", porcentaje, minutos, segundos_restantes);
    fflush(stdout);
}

// Función para inicializar todas las nuevas grids al comienzo de una simulación.

void initialize_coarse_grid() {

    memset(W_grid, 0, sizeof(W_grid));
    next_k_value = 1;
    memset(Nk_counts, 0, sizeof(Nk_counts));

}

void add_particle_to_coarse_grid(int p_index) {

    double px = particles_list[p_index][0];
    double py = particles_list[p_index][1];

    int cell_x = (int)floor(px/CELL_SIZE);
    int cell_y = (int)floor(py/CELL_SIZE);

    if (cell_x < 0 || cell_x >= GRID_DIM_X || cell_y < 0 || cell_y >= GRID_DIM_Y) return;

    int k = W_grid[cell_x][cell_y];

    if (k == 0) {

        k = next_k_value;
        W_grid[cell_x][cell_y] = k;
        next_k_value++;
    }

    int count = Nk_counts[k];
    
    // Reactivar la verificación de desbordamiento

    if (count >= BLOCK_SIZE) {

        fprintf(stderr, "\nError: La celda de la grilla gruesa (%d, %d) está llena. Aumente L_PARAM.\n", cell_x, cell_y);
        exit(1); // Detiene el programa si una celda se desborda.
    }
    
    int start_of_block = BLOCK_SIZE * (k - 1);
    F_indices[start_of_block + count] = p_index;
    Nk_counts[k]++;
}

void inicializar_grids_optimizacion() {

    // Calcular la Vicinity Grid (Ψ) una sola vez

    for (int i = 0; i < 2 * D_MAX + 1; i++) {

        for (int j = 0; j < 2 * D_MAX + 1; j++) {

            double dist = sqrt((i - D_MAX)*(i - D_MAX) + (j - D_MAX)*(j - D_MAX));
            Vicinity_grid[i][j] = (int)round(dist);
        }
    }

    // Inicializar Grids Y y Ω

    for (int i = 0; i < (int)SYSTEM_WIDTH; i++) {

        for (int j = 0; j < Lm; j++) {

            Omega_grid[i][j] = D_MAX;
            Y_grid[i][j] = 0; // 0 significa vacío
        }
    }
}

void add_and_update_all_grids(int p_index) {

    // Añadir a la grid gruesa para predicción de colisión

    add_particle_to_coarse_grid(p_index);

    // Actualizar Y_grid y Omega_grid para la optimización de trayectoria

    double px = particles_list[p_index][0];
    double py = particles_list[p_index][1];
    
    int px_grid = (int)round(px);
    int py_grid = (int)round(py);

    // Aplicar periodicidad a la coordenada X

    px_grid = (px_grid % (int)SYSTEM_WIDTH + (int)SYSTEM_WIDTH)%(int)SYSTEM_WIDTH;

    if (py_grid < 0 || py_grid >= Lm) return;
    
    Y_grid[px_grid][py_grid] = p_index + 1; // Guardar índice + 1

    // Actualizar Omega_grid alrededor de la nueva partícula

    for (int i = -D_MAX; i <= D_MAX; i++) {

        for (int j = -D_MAX; j <= D_MAX; j++) {

            int gx = px_grid + i;
            int gy = py_grid + j;

            // Aplicar periodicidad a la coordenada x de la grid

            gx = (gx % (int)SYSTEM_WIDTH + (int)SYSTEM_WIDTH)%(int)SYSTEM_WIDTH;

            if (gy >= 0 && gy < Lm) {
                int d = Vicinity_grid[i + D_MAX][j + D_MAX];
                if (d < Omega_grid[gx][gy]) Omega_grid[gx][gy] = d;
            }
        }
    }
}

double solve_quadratic(double a, double b, double c) {

    double discriminant = b*b - 4.0*a*c;
    if (discriminant < 0) return -1.0;

    double sqrt_discriminant = sqrt(discriminant);
    double root1 = (-b + sqrt_discriminant)/(2.0 * a);
    double root2 = (-b - sqrt_discriminant)/(2.0 * a);

    double epsilon = 1e-9;
    if (root1 > epsilon && root2 > epsilon) return fmin(root1, root2);
    if (root1 > epsilon) return root1;
    if (root2 > epsilon) return root2;
    return -1.0;
}

double find_collision_distance(double walker_x, double walker_y, double move_theta, double L_max) {

    double min_L_hit = L_max + 1.0;
    
    int cell_x = (int)floor(walker_x/CELL_SIZE);
    int cell_y = (int)floor(walker_y/CELL_SIZE);
    
    double cos_theta = cos(move_theta);
    double sin_theta = sin(move_theta);

    for (int dx = -1; dx <= 1; dx++) {

        for (int dy = -1; dy <= 1; dy++) {

            int check_cell_x = (cell_x + dx + GRID_DIM_X)%GRID_DIM_X; // Periodicidad
            int check_cell_y = cell_y + dy;

            if (check_cell_y >= 0 && check_cell_y < GRID_DIM_Y) {

                int k = W_grid[check_cell_x][check_cell_y];
                if (k == 0) continue;

                int num_particles_in_cell = Nk_counts[k];
                int start_of_block = BLOCK_SIZE*(k - 1);

                for (int i = 0; i < num_particles_in_cell; i++) {

                    int p_index = F_indices[start_of_block + i];
                    double p_x = particles_list[p_index][0];
                    double p_y = particles_list[p_index][1];

                    // Calcular la distancia periódica más corta en x
                    double dist_x = walker_x - p_x;
                    if (dist_x > SYSTEM_WIDTH/2.0) dist_x -= SYSTEM_WIDTH;
                    if (dist_x < -SYSTEM_WIDTH/2.0) dist_x += SYSTEM_WIDTH;

                    double dist_y = walker_y - p_y;

                    // Optimización: salto rápido si la partícula está demasiado lejos

                    double max_safe_dist = L_max + diameter;
                    if (fabs(dist_x) > max_safe_dist || fabs(dist_y) > max_safe_dist) continue;

                    // Coeficientes para la ecuación cuadrática

                    double B = 2.0*(cos_theta*dist_x + sin_theta*dist_y);
                    double C = dist_x*dist_x + dist_y*dist_y - collision_dist_sq;
                    double L_hit = solve_quadratic(1.0, B, C);

                    if (L_hit > 1e-9 && L_hit <= L_max && L_hit < min_L_hit) min_L_hit = L_hit;
                }
            }
        }
    }

    return (min_L_hit <= L_max) ? min_L_hit : -1.0;
}

void save_particles(const char* filename, int count) {

    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("No se pudo abrir el archivo");
        return;
    }
    for (int i = 0; i < count; i++) {
        fprintf(file, "%.4f %.4f\n", particles_list[i][0], particles_list[i][1]);
    }
    fclose(file);
}

void get_interface_profile(int p_count, double* out_heights) {

    for (int i = 0; i < L; i++) out_heights[i] = 0.0;

    for (int i = 0; i < p_count; i++) {

        double px = particles_list[i][0];
        double py = particles_list[i][1];

        int bin_index = (int)floor(px / diameter);

        if (bin_index >= 0 && bin_index < L) {

            if (py > out_heights[bin_index]) {

                out_heights[bin_index] = py;
            }
        }
    }
}

// Espesor RMS

double get_rms_thickness(int current_particle_count) {
    
    // Si solo tenemos el sustrato inicial, el espesor del depósito es 0

    if (current_particle_count <= L) return 0.0;

    long double sum_y = 0.0;
    long double sum_y_sq = 0.0;
    int count_deposited = 0;

    // Recorremos SOLO las partículas depositadas
    // (Desde L en adelante)

    for (int i = L; i < current_particle_count; i++) {
        
        long double y = particles_list[i][1];
        
        sum_y += y;
        sum_y_sq += y * y;
        count_deposited++;
    }

    if (count_deposited == 0) return 0.0;

    long double N_rms = (long double)count_deposited;
    long double mean_y = sum_y / N_rms;
    long double mean_y_sq = sum_y_sq / N_rms;

    // Varianza = <y^2> - <y>^2

    long double variance = mean_y_sq - (mean_y * mean_y);

    // RMS Thickness T = raiz(Varianza)
    // Protección contra errores de punto flotante negativos muy pequeños

    if (variance < 0.0) return 0.0;
    
    return sqrt((double)variance);
}

// Calculo de la dimensión fractal

// Función auxiliar: Comprueba si hay ALGO en la caja

int any_box(int Y_grid[][Lm], int r_start, int c_start, int size, int limit_x, int limit_y, double particles[][2], double diam) {
    
    int r_end = r_start + size;
    int c_end = c_start + size;
    
    if (c_end > limit_y) c_end = limit_y;
    
    double radius = diam/2.0;
    double box_x_min = (double)r_start;
    double box_x_max = (double)(r_start + size);
    double box_y_min = (double)c_start;
    double box_y_max = (double)(c_start + size);
    
    for (int i = r_start; i < r_end; i++) {
        
        // Aplicar periodicidad correctamente

        int grid_i = i % limit_x;
        if (grid_i < 0) grid_i += limit_x;
        
        for (int j = c_start; j < c_end; j++) {
            
            int p_id = Y_grid[grid_i][j];
            
            if (p_id > 0) {

                int p_index = p_id - 1;
                double p_real_x = particles[p_index][0];
                double p_real_y = particles[p_index][1];
                
                // Verificación en Y (simple, no periódico)

                int fits_y = (p_real_y - radius >= box_y_min) && 
                            (p_real_y + radius <= box_y_max);
                
                if (!fits_y) continue;
                
                // Verificación en X con periodicidad
                // Calcular distancia mínima considerando wrapping

                double dx_min = fabs(p_real_x - box_x_min);
                double dx_max = fabs(p_real_x - box_x_max);
                
                if (dx_min > limit_x/2.0) dx_min = limit_x - dx_min;
                if (dx_max > limit_x/2.0) dx_max = limit_x - dx_max;
                
                // La partícula cabe si está suficientemente dentro

                int fits_x = (dx_min >= radius) && (dx_max >= radius);
                
                if (fits_x && fits_y) {
                    return 1;
                }
            }
        }
    }
    
    return 0;
}

// Función principal

int f_dimension(int grid[][Lm], int max_height, long double results[30][2], int size_L, double particles[][2], double diam) {
    
    int step_index = 0;
    
    // Calcular límites reales

    int actual_width = SYSTEM_WIDTH;  // Ancho real de la grid
    int limit_y = (int)ceil(max_height) + 2;  // Altura + margen
    if (limit_y > Lm) limit_y = Lm;
    
    // Empezar con tamaño razonable (potencia de 2 menor a max_height)

    int size = 2;

    while (size*2 < limit_y && size*2 < actual_width) {

        size *= 2;
    }
    
    int min_size = 2*diameter;
    
    while (size >= min_size && step_index < 30) {
        
        long long count = 0;
        
        // Barrido correcto

        for (int i = 0; i < actual_width; i += size) {

            for (int j = 1; j < limit_y; j += size) {
                
                if (any_box(Y_grid, i, j, size, actual_width, limit_y, particles_list, diam)) {

                    count++;
                }
            }
        }
        
        results[step_index][0] = (long double)size;
        results[step_index][1] = (long double)count;
        
        size /= 2;
        step_index++;
    }
    
    return step_index;
}

// ############################################ Inicio de la simulación #########################################################

int main(int argc, char *argv[]) {

    job_id = 0; // ID por defecto

    if (argc > 1) {
        job_id = atoi(argv[1]); // Capturar el ID desde la línea de comandos
    }

    // Garantizar que incluso si dos trabajos empiezan en el mismo segundo, la semilla será diferente.

    xorshift_seed(time(NULL) ^ getpid() ^ job_id); 
    
    time_t tiempo_inicio = time(NULL);

    collision_dist_sq = diameter * diameter;

    // Inicialización de vectores y matrices

    memset(M1, 0, sizeof(M1));
    memset(M2, 0, sizeof(M2));
    memset(M3, 0, sizeof(M3));
    memset(M4, 0, sizeof(M4));
    memset(TIME, 0, sizeof(TIME));
    memset(N_PARTICLES, 0, sizeof(N_PARTICLES));

    remove("all_trees_raw.dat");

    memset(H_sum_binned, 0, sizeof(H_sum_binned));
    memset(Count_binned, 0, sizeof(Count_binned));

    memset(total_rms_per_layer, 0, sizeof(total_rms_per_layer));
    
    inicializar_grids_optimizacion();
    init_tree_analysis();  // Inicializar arrays de promedio

    for (unsigned int m = 0; m < N; m++) {

        mostrar_barra_progreso(m + 1, N, "Procesando Muestra", tiempo_inicio);

        initialize_coarse_grid();

        // Las grids Y y Omega se resetean lógicamente al crear el sustrato

        for (int i = 0; i < (int)SYSTEM_WIDTH; i++) {

            for (int j = 0; j < Lm; j++) {

                Omega_grid[i][j] = D_MAX;
                Y_grid[i][j] = 0;
            }
        }
        
        particle_count = 0;
        double sample_max_h = diameter; 

        // Crear el sustrato inicial

        for (int i = 0; i < L; i++) {

            particles_list[particle_count][0] = (double)i*diameter;
            particles_list[particle_count][1] = 0.0;
            add_and_update_all_grids(particle_count);
            particle_count++;
        }

        // Deposición

        for (int t = 0; t < T; t++) {

            int active_h_count = 0;
            int start_saving_index = L - LAST_PARTICLES_TO_SAVE;

            for (int i = 0; i < L; i++) {

                double hmax = sample_max_h + 100;

                double walker_x = rand_double()*SYSTEM_WIDTH;
                double walker_y = sample_max_h + hmax; // Se lanza desde la altura diameter + 100

                int stuck = 0;

                while (!stuck) {

                    double current_step;
                    
                    int walker_ix = (int)round(walker_x);
                    int walker_iy = (int)round(walker_y);

                    // Aplicar periodicidad para la consulta en la grid

                    walker_ix = (walker_ix%(int)SYSTEM_WIDTH + (int)SYSTEM_WIDTH)%(int)SYSTEM_WIDTH;
                    
                    int d_wc = D_MAX;

                    if (walker_iy >= 0 && walker_iy < Lm) d_wc = Omega_grid[walker_ix][walker_iy];

                    if (d_wc <= COLLISION_THRESHOLD) current_step = L_MIN;
                    else current_step = d_wc - COLLISION_THRESHOLD
                    double move_theta = rand_double()*2.0*M_PI;

                    if (current_step == L_MIN) {

                        double L_hit = find_collision_distance(walker_x, walker_y, move_theta, L_MIN);

                        if (L_hit > 0) { // Colisión

                            walker_x += L_hit*cos(move_theta);
                            walker_y += L_hit*sin(move_theta);
                            stuck = 1;

                        } else { // No hay colisión, dar paso completo

                            walker_x += current_step*cos(move_theta);
                            walker_y += current_step*sin(move_theta);
                        }

                    } else { // Salto grande y seguro

                        walker_x += current_step*cos(move_theta);
                        walker_y += current_step*sin(move_theta);
                    }
                    
                    // Aplicar condiciones de contorno periódicas al caminante

                    walker_x = fmod(walker_x, SYSTEM_WIDTH);
                    if (walker_x < 0) walker_x += SYSTEM_WIDTH;

                    if (stuck) {

                        particles_list[particle_count][0] = walker_x;
                        particles_list[particle_count][1] = walker_y;
                        add_and_update_all_grids(particle_count);

                        if (walker_y > sample_max_h) sample_max_h = walker_y;
                        particle_count++;

                        if (i >= start_saving_index) {

                            if (active_h_count < LAST_PARTICLES_TO_SAVE) {
                                active_h_vector[active_h_count] = walker_y; 
                                active_h_count++;
                            }
                        }

                    } else {

                        // Condición de escape

                        if (walker_y > 5*sample_max_h || walker_y < -diameter) {

                            // Partícula perdida, relanzar

                            walker_x = rand_double()*SYSTEM_WIDTH;
                            walker_y = sample_max_h + hmax;
                        }
                    }

                } // Fin while(!stuck)

            } // Fin bucle deposición por capa

            get_interface_profile(particle_count, interface_heights);

            N_PARTICLES[t] += (long double)particle_count;

            long double hm1 = 0.0, hm2 = 0.0, hm3 = 0.0, hm4 = 0.0;

            for (int k = 0; k < L; k++) {

                long double hi = (long double)interface_heights[k];

                hm1 += hi;
                hm2 += hi*hi;
                hm3 += hi*hi*hi;
                hm4 += hi*hi*hi*hi;
            }

            M1[t] += hm1/(long double)L;
            M2[t] += hm2/(long double)L;
            M3[t] += hm3/(long double)L;
            M4[t] += hm4/(long double)L;

            long double h_act1 = 0.0, h_act2 = 0.0, h_act3 = 0.0, h_act4 = 0.0;

            for (int h_idx = 0; h_idx < active_h_count; h_idx++) {

                long double hi_act = (long double)active_h_vector[h_idx];

                h_act1 += hi_act;
                h_act2 += hi_act*hi_act;
                h_act3 += hi_act*hi_act*hi_act;
                h_act4 += hi_act*hi_act*hi_act*hi_act;
            }


            M1_active[t] += h_act1/(long double)active_h_count;
            M2_active[t] += h_act2/(long double)active_h_count;
            M3_active[t] += h_act3/(long double)active_h_count;
            M4_active[t] += h_act4/(long double)active_h_count;

            double current_rms = get_rms_thickness(particle_count);
            total_rms_per_layer[t] += current_rms;


            TIME[t] = t + 1;

            int current_t = t + 1;
            int snapshot_found = -1;

            // Verificamos si el tiempo actual coincide con alguno de nuestros objetivos

            for(int k=0; k < NUM_SNAPSHOTS; k++){

                if(current_t == snapshot_times[k]){

                    snapshot_found = k;
                    break;
                }
            }

            // Si es un tiempo objetivo, ejecutamos el análisis

            if (snapshot_found != -1) {

                int snapshot_idx = snapshot_found; 

                analyze_trees_bfs(particle_count, diameter, m, current_t, snapshot_found);

                memset(results, 0, sizeof(results));

                int pasos = f_dimension(Y_grid, (int)sample_max_h, results, L, particles_list, diameter);

                for (int k = 0; k < pasos; k++) {

                    total_results[snapshot_idx][k][0] += results[k][0]; // Guardar en su "slot" de tiempo
                    total_results[snapshot_idx][k][1] += results[k][1];
                }
            }

        } // Fin bucle de tiempo

        if (sample_max_h >= Lm) {

            printf("\n\nADVERTENCIA: La altura máximo (%f) ha alcanzado el borde de la grilla (Lm = %d).\n", sample_max_h, Lm);
            printf("Simulación detenida en la muestra %d de %d.\n", m + 1, N);
            printf("Se guardarán los datos recopilados hasta este punto.\n");
            break;
        }

    } // Fin bucle de muestras

    char rms_filename[60];

    sprintf(rms_filename, "rms_growth_avg_L=%d_run=%d.dat", L, job_id);
    FILE* rms_file = fopen(rms_filename, "w");

    if (rms_file == NULL) {
        perror("Error creando archivo RMS");
    } else {
        fprintf(rms_file, "# N_Depositadas\tT_rms_avg\n");

        for (int i = 0; i < T; i++) {
            
            // Calculamos N teórico promedio por capa (L partículas por capa)
            // Esto es más limpio que usar particle_count que se reinicia.

            long long N_acumulado = (long long)(i + 1)*L;
            
            // Promedio sobre las N muestras

            double avg_T = (double)(total_rms_per_layer[i]/(long double)N);
            
            fprintf(rms_file, "%lld\t%.6f\n", N_acumulado, avg_T);
        }

        fclose(rms_file);
    }

    printf("\nGuardando análisis de interfaz...\n");
    
    char file_w2[60], file_sk[60], file_kt[60], file_mean[60];
    char file_w2_active[60], file_sk_active[60], file_kt_active[60], file_mean_active[60];

    sprintf(file_w2, "roughness_opt_L=%d_run=%d.dat", L, job_id);
    sprintf(file_sk, "skewness_opt_L=%d_run=%d.dat", L, job_id);
    sprintf(file_kt, "kurtosis_opt_L=%d_run=%d.dat", L, job_id);
    sprintf(file_mean, "mean_height_opt_L=%d_run=%d.dat", L, job_id);

    sprintf(file_w2_active, "roughness_active_opt_L=%d_run=%d.dat", L, job_id);
    sprintf(file_sk_active, "skewness_active_opt_L=%d_run=%d.dat", L, job_id);
    sprintf(file_kt_active, "kurtosis_active_opt_L=%d_run=%d.dat", L, job_id);
    sprintf(file_mean_active, "mean_height_active_opt_L=%d_run=%d.dat", L, job_id);

    FILE *filevar = fopen(file_w2, "w"), *filesk = fopen(file_sk, "w"), *filekt = fopen(file_kt, "w"), 
    *filemean_h = fopen(file_mean, "w");

    FILE *filevar_active = fopen(file_w2_active, "w"), *filesk_active = fopen(file_sk_active, "w"), *filekt_active = fopen(file_kt_active, "w"), 
    *filemean_h_active = fopen(file_mean_active, "w");

    for (int tm = 0; tm < T; tm++) {

        M1[tm] /= (long double)N;
        M2[tm] /= (long double)N;
        M3[tm] /= (long double)N;
        M4[tm] /= (long double)N;

        long double varC = M2[tm] - (M1[tm] * M1[tm]);
        long double k3C = M3[tm] - 3*M1[tm]*M2[tm] + 2*(M1[tm]*M1[tm]*M1[tm]);
        long double k4C = M4[tm] - 4*M1[tm]*M3[tm] - 3*(M2[tm]*M2[tm]) + 12*(M1[tm]*M1[tm])*M2[tm] - 6*(M1[tm]*M1[tm]*M1[tm]*M1[tm]);

        long double skewnessC = (varC > 1e-9) ? k3C/(varC * sqrt(varC)) : 0.0;
        long double kurtosisC = (varC > 1e-9) ? k4C/(varC * varC) : 0.0;

        fprintf(filevar, "%d %.9Lf\n", TIME[tm], varC);
        fprintf(filesk, "%d %.9Lf\n", TIME[tm], skewnessC);
        fprintf(filekt, "%d %.9Lf\n", TIME[tm], kurtosisC);

        fprintf(filemean_h, "%d %.9Lf\n", TIME[tm], M1[tm]);

        M1_active[tm] /= (long double)N;
        M2_active[tm] /= (long double)N;
        M3_active[tm] /= (long double)N;
        M4_active[tm] /= (long double)N;

        long double varC_active = M2_active[tm] - (M1_active[tm] * M1_active[tm]);
        long double k3C_active = M3_active[tm] - 3*M1_active[tm]*M2_active[tm] + 2*(M1_active[tm]*M1_active[tm]*M1_active[tm]);
        long double k4C_active = M4_active[tm] - 4*M1_active[tm]*M3_active[tm] - 3*(M2_active[tm]*M2_active[tm]) 
                            + 12*(M1_active[tm]*M1_active[tm])*M2_active[tm] 
                            - 6*(M1_active[tm]*M1_active[tm]*M1_active[tm]*M1_active[tm]);

        long double skewnessC_active = (varC_active > 1e-9) ? k3C_active/(varC_active * sqrt(varC_active)) : 0.0;
        long double kurtosisC_active = (varC_active > 1e-9) ? k4C_active/(varC_active * varC_active) : 0.0;

        fprintf(filevar_active, "%d %.9Lf\n", TIME[tm], varC_active);
        fprintf(filesk_active, "%d %.9Lf\n", TIME[tm], skewnessC_active);
        fprintf(filekt_active, "%d %.9Lf\n", TIME[tm], kurtosisC_active);

        fprintf(filemean_h_active, "%d %.9Lf\n", TIME[tm], M1_active[tm]);

    }

    fclose(filevar);
    fclose(filesk);
    fclose(filekt);
    fclose(filemean_h);
    fclose(filevar_active);
    fclose(filesk_active);
    fclose(filekt_active);
    fclose(filemean_h_active);

    for (int s = 0; s < NUM_SNAPSHOTS; s++) {

        char archivo_frac[100];
        int t_val = snapshot_times[s];

        sprintf(archivo_frac, "fractal_dim_t=%d_L=%d_run=%d.dat", t_val, L, job_id);
        
        FILE *f_frac = fopen(archivo_frac, "w");
        fprintf(f_frac, "# Size(s)  Count(N) (Promedio sobre %d muestras)\n", N);

        for (int k = 0; k < 30; k++) {

            long double avg_size = total_results[s][k][0]/(long double)N; // Usar índice [s]
            long double avg_count = total_results[s][k][1]/(long double)N;

            if (avg_size == 0) break; // Salir si no hay más datos válidos

            fprintf(f_frac, "%.9Lf %.9Lf\n", avg_size, avg_count);
        }

        fclose(f_frac);
    }

    save_averaged_tree_data();

    save_particles("particulas_final.dat", particle_count);
    
    printf("Simulación finalizada. Coordenadas y análisis guardados.\n");

    return 0;
}

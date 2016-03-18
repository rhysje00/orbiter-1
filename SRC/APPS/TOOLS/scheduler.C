// scheduler.C
// 
// Anton Betten
// January 16, 2016
//
// 
//

#include "orbiter.h"
#include "discreta.h"

typedef class job_table job_table;

class job_table {

public:

	job_table();
	~job_table();

	INT f_task_assigned;
	INT task;
	INT job;
	BYTE target_fname[1000];
	BYTE command[1000];
	const BYTE *target_file_mask;
	const BYTE *command_mask;
};


// global data:

INT t0; // the system time when the program started

void do_scheduling(INT N, INT J, const BYTE *target_file_mask, 
	INT f_command_mask, const BYTE *command_mask, 
	INT *excluded_cases, INT nb_excluded_cases,
	INT f_reload, 
	INT verbose_level);
INT find_free_job(job_table *JT, INT J);
void assign_task(job_table *JT, INT t, INT j, INT verbose_level);


int main(int argc, char **argv)
{
	INT i;
	t0 = os_ticks();
	INT verbose_level = 0;
	INT f_target_file_mask = FALSE;	
	const BYTE *target_file_mask = NULL;
	INT f_command_mask = FALSE;	
	const BYTE *command_mask = NULL;
	INT f_N = FALSE;
	INT N = 0;
	INT f_J = FALSE;
	INT J = 1;
	INT nb_excluded_cases = 0;
	INT excluded_cases[1000];
	INT f_reload = FALSE;

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-v") == 0) {
			verbose_level = atoi(argv[++i]);
			cout << "-v " << verbose_level << endl;
			}
		else if (strcmp(argv[i], "-target_file_mask") == 0) {
			f_target_file_mask = TRUE;
			target_file_mask = argv[++i];
			cout << "-target_file_mask " << target_file_mask << endl;
			}
		else if (strcmp(argv[i], "-command_mask") == 0) {
			f_command_mask = TRUE;
			command_mask = argv[++i];
			cout << "-command_mask " << command_mask << endl;
			}
		else if (strcmp(argv[i], "-N") == 0) {
			f_N = TRUE;
			N = atoi(argv[++i]);
			cout << "-N " << N << endl;
			}
		else if (strcmp(argv[i], "-J") == 0) {
			f_J = TRUE;
			J = atoi(argv[++i]);
			cout << "-J " << J << endl;
			}
		else if (strcmp(argv[i], "-exclude") == 0) {
			excluded_cases[nb_excluded_cases] = atoi(argv[++i]);
			cout << "-exclude " << excluded_cases[nb_excluded_cases] << endl;
			nb_excluded_cases++;
			}
		else if (strcmp(argv[i], "-reload") == 0) {
			f_reload = TRUE;
			cout << "-reload " << endl;
			}
		}

	if (!f_N) {
		cout << "please use option -N <N>" << endl;
		exit(1);
		}
	if (!f_J) {
		cout << "please use option -J <J>" << endl;
		exit(1);
		}
	if (!f_target_file_mask) {
		cout << "please use option -target_file_mask <target_file_mask>" << endl;
		exit(1);
		}

	INT_vec_heapsort(excluded_cases, nb_excluded_cases);
	cout << "There are " << nb_excluded_cases << " excluded cases: ";
	INT_vec_print(cout, excluded_cases, nb_excluded_cases);
	cout << endl;
	

	do_scheduling(N, J, target_file_mask, 
		f_command_mask, command_mask, 
		excluded_cases, nb_excluded_cases,
		f_reload, 
		verbose_level);

	cout << "scheduler.out is done" << endl;
	the_end(t0);
	//the_end_quietly(t0);

}

void do_scheduling(INT N, INT J, const BYTE *target_file_mask, 
	INT f_command_mask, const BYTE *command_mask, 
	INT *excluded_cases, INT nb_excluded_cases,
	INT f_reload, 
	INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	INT i;

	if (f_v) {
		cout << "do_scheduling" << endl;
		cout << "N = " << N << endl;
		cout << "J = " << J << endl;
		cout << "target_file_mask = " << target_file_mask << endl;
		cout << "f_command_mask = " << f_command_mask << endl;
		if (f_command_mask) {
			cout << "command_mask = " << command_mask << endl;
			}
		cout << "f_reload = " << f_reload << endl;
		}

	INT *task_completed;
	INT nb_tasks_completed = 0;
	BYTE target_fname[1000];

	task_completed = NEW_INT(N);

	for (i = 0; i < N; i++) {
		sprintf(target_fname, target_file_mask, i);
		if (file_size(target_fname) > 0) {
			task_completed[i] = 1;
			nb_tasks_completed++;
			}
		else {
			task_completed[i] = 0;
			}
		}
	
	cout << "number of completed tasks = " << nb_tasks_completed << " / " << N << endl;

	cout << "there are " << N - nb_tasks_completed << " open tasks" << endl;
	for (i = 0; i < N; i++) {
		if (task_completed[i] == 0) {
			cout << i << " ";
			}
		}
	cout << endl;


	INT idx;
	
	if (f_command_mask) {

		job_table *JT;

		JT = new job_table[J];
		for (i = 0; i < J; i++) {
			JT[i].job = i;
			JT[i].f_task_assigned = FALSE;
			JT[i].target_file_mask = target_file_mask;
			JT[i].command_mask = command_mask;
			}

		while (TRUE) {

			INT t, j;
	
			for (t = 0; t < N; t++) {
				if (task_completed[t] == 0 && !INT_vec_search(excluded_cases, nb_excluded_cases, t, idx)) {
					j = find_free_job(JT, J);
					if (j == -1) {
						break;
						}
					assign_task(JT, t, j, verbose_level);
					task_completed[t] = 2;
					}
				}

			if (!f_reload) {
				break;
				}

			system("sleep 5");


			for (i = 0; i < J; i++) {

				if (!JT[i].f_task_assigned) {
					continue;
					}
				t = JT[i].task;
				sprintf(target_fname, target_file_mask, t);
				if (file_size(target_fname) > 0) {
					cout << "task completed: " << t << endl;
					task_completed[t] = 1;
					JT[i].f_task_assigned = FALSE;
					nb_tasks_completed++;
					}
				}
	
			cout << "number of completed tasks = " << nb_tasks_completed << " / " << N << endl;

			cout << "there are " << N - nb_tasks_completed << " open tasks" << endl;

			if (nb_tasks_completed == N) {
				break;
				}

			if (N - nb_tasks_completed < 100) {
				for (i = 0; i < N; i++) {
					if (task_completed[i] == 0) {
						cout << i << " ";
						}
					}
				cout << endl;
				}
			else {
				cout << "too many to print" << endl;
				}




			}
		}
}

INT find_free_job(job_table *JT, INT J)
{
	INT i;

	for (i = 0; i < J; i++) {
		if (!JT[i].f_task_assigned) {
			return i;
			}
		}
	return -1;
}

void assign_task(job_table *JT, INT t, INT j, INT verbose_level)
{
	INT f_v = (verbose_level >= 1);
	BYTE str[1000];

	if (f_v) {
		cout << "assign_task assigning task " << t << " to job " << j << endl;
		}
	sprintf(JT[j].target_fname, JT[j].target_file_mask, t);
	sprintf(JT[j].command, JT[j].command_mask, t, t, t);
	sprintf(str, "  >log_%ld &", t);
	strcat(JT[j].command, str);
	JT[j].task = t;
	JT[j].f_task_assigned = TRUE;
	cout << "assign task: " << JT[j].command << endl;
	system(JT[j].command);
}

job_table::job_table()
{
}

job_table::~job_table()
{
}



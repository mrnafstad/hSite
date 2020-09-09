<template>
	<div>
		<div v-if="loaded">
			<h3>Todos</h3>
			<div class="legend">
				<span>Click to mark as complete</span>
				<span>
					<span class="incomplete-box"></span> = Incomplete
				</span>
				<span>
					<span class="complete-box"></span> = Complete
				</span>
			</div>
			<AddTodo />
			<FilterTodos />
			<div class="todos">
				<div @click="onDblClick(todo)" v-for="(todo, idx) in allTodos" :key="idx" class="todo" v-bind:class="{'is-complete':todo.completed}">
					{{ todo.title }}
					<i @click="deleteTodo(todo.id)" class="fas fa-trash-alt"></i>
				</div>
			</div>
		</div>
		<div v-if="!loaded"> Loading todos... </div>
	</div>
</template>

<script>
import { mapGetters, mapActions } from 'vuex'
import FilterTodos from '@/components/FilterTodos.vue'
import AddTodo from '@/components/AddTodo.vue'

export default {
	data() {
		return {
			name: "Todos",
			loading: true,
			loaded: false
		}
	},
	components: {
		FilterTodos,
		AddTodo
	},
	methods: {
		...mapActions(['fetchTodos', "deleteTodo", "updateTodo"]),
		onDblClick(todo) {
			const updTodo = {
				id: todo.id,
				title: todo.title,
				completed: !todo.completed
			}

			this.updateTodo(updTodo)
		}
	},
	computed: mapGetters(['allTodos']),
	created() {
		this.fetchTodos().then(
			this.loaded = true)
	}
}
</script>

<style scoped>
	.todos {
		display: grid;
		grid-template-columns: repeat(3, 1fr);
		grid-gap: 1rem;
		padding: 1rem;
	}

	.todo {
		border: 1px solid #ccc;
		background: #41b883;
		padding: 1rem;
		border-radius: 5px;
		text-align: center;
		positon: relative;
		cursor: pointer;
	}

	i {
		/*position: absolute;*/
		bottom: 10px;
		right: 10px;
		color: #fff;
		cursor: pointer;
	}

	.legend {
		display: flex;
		justify-content: space-around;
		margin-bottom: 1rem;
	}

	.complete-box {
		display: inline-block;
		width: 15px;
		height: 10px;
		background: #35495e;
	}

	.incomplete-box {
		display: inline-block;
		width: 15px;
		height: 10px;
		background: #41b883;		
	}

	@media (max-width: 500px) {
		.todos {
			grid-template-columns: 1fr;
		}
	}

	.is-complete {
		background: #35405e;
		color: #fff;
	}
</style>

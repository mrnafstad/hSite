<template>
	<div>
		<h3>Todos</h3>
		<div class="legend">
			<span>Double click to mark as complete</span>
			<span>
				<span class="incomplete-box"></span> = Incomplete
			</span>
			<span>
				<span class="complete-box"></span> = Complete
			</span>
		</div>
		<div class="todos">
			<div @click="onDblClick(todo)" v-for="todo in allTodos" :key="todo.id" class="todo" v-bind:class="{'is-complete':todo.completed}">
				{{ todo.title }}
				<i id="trash" @click="deleteTodo(todo.id)" class="fas fa-trash-alt"></i>
			</div>
		</div>
	</div>
</template>

<script>
import { mapGetters, mapActions } from 'vuex'

export default {
	name: "Todos",
	methods: {
		...mapActions(['fetchTodos', "deleteTodo", "updateTodo"]),
		onDblClick(todo) {
			console.log(todo.id)
			const updTodo = {
				id: todo.id,
				title: todo.title,
				completed: !todo.completed
			}
			console.log(updTodo)

			this.updateTodo(updTodo)
		}
	},
	computed: mapGetters(['allTodos']),
	created() {
		this.fetchTodos();
	}
}
</script>

<style scoped>
	.todos {
		display: grid;
		grid-template-columns: repeat(3, 1fr);
		grid-gap: 1rem;
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

	#trash {
		position: absolute;
		bottom: 10rem;
		right: 10rem;
		color: black;
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

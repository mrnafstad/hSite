import axios from 'axios'
import db from '../../firebaseConfig.js'
//const db = require('../../firebaseConfig.js')

const state = {
	isFetched: false,
	todos: []
}

const getters = {
	allTodos: (state) => state.todos
}

const actions = {
	async fetchTodos({ commit }) {
		if (!state.isFetched) {
			await db.todos.get().then(
				function(querySnapshot) {
					querySnapshot.forEach(function(doc) {
						//console.log(doc.id, "=>", doc.data())
						commit('newTodo', {
							id: doc.id,
							title: doc.data().title,
							completed: doc.data().completed
						})
						state.isFetched = true
					})
					//commit('setTodos', querySnapshot)
					console.log("Todos retrieved from database: ", state.isFetched)
				})
			//const response = await axios.get('https://jsonplaceholder.typicode.com/todos')
			//commit('setTodos', response.data)
		}
	},
	async addTodo({ commit }, title) {
		//const response = await axios.post('https://jsonplaceholder.typicode.com/todos', {title, completed: false})
		await db.todos.add({
			title,
			completed: false
		})
		.then(function(docRef) {
			commit('newTodo', docRef)
			//console.log("Document successfully written =>", docRef.id)
		})
		.catch(function(error) {
			console.error("Error writing document: ", error)
		})
	},
	async deleteTodo({ commit }, id) {
		//await axios.delete(`https://jsonplaceholder.typicode.com/todos/${id}`)

		await db.todos.doc(id).delete().then(function() {
			console.log("Todo deleted!")
			commit('removeTodo', id)
		}).catch(function(error) {
			console.error("Error removing todo: ", error)
		})
	},
	async filterTodos({ commit }, e) {
		const limit = parseInt(e.target.value)

		const response = await axios.get(`https://jsonplaceholder.typicode.com/todos?_limit=${limit}`)

		commit('setTodos', response.data)
	},
	async updateTodo({ commit }, updTodo) {
		//const response = await axios.put(`https://jsonplaceholder.typicode.com/todos/${updTodo.id}`, updTodo)
		await db.todos.doc(updTodo.id).update({
			completed: updTodo.completed
		})
		.then(function() {
			commit('updateTodo', updTodo)
			console.log("Updated ID: ", updTodo.id)
		})
		.catch(function(err) {
			console.log(updTodo)
			console.error("Error updating element: ", err)
		})
	}
}

const mutations = {
	setTodos: (state, todos) => (state.todos = todos),
	newTodo: (state, todo) => state.todos.unshift(todo),
	removeTodo: (state, id) => state.todos = state.todos.filter(todo => todo.id !== id),
	updateTodo: (state, updTodo) => {
		const index = state.todos.findIndex(todo => todo.id === updTodo.id)
		if(index !== -1) {
			state.todos.splice(index, 1, updTodo)
		}
	}
}

export default {
	state,
	getters,
	actions,
	mutations
}

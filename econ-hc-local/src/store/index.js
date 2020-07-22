import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
	state: {
		bruker: {},
		andreBrukere: [],
		loggin: false
	},
	mutations: {
		setPerson (state, person) {
			if (!state.login) {
				state.bruker = person
				state.loggin = true
				console.log(state.bruker)
			}
		},
		setEksterneBrukere(state, brukere) {
			state.andreBrukere.push(brukere)
			//console.log(state.andreBrukere)
		}
	},
	getters: {
		bruker: state => {
			if (state.loggin) {
				return state.bruker
			}
		},
		login: state => {
			return state.loggin
		},
		andre: state => {
			return state.andreBrukere
		}
	},
	actions: {
	},
	modules: {
	}
})

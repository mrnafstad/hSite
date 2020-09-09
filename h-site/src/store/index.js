import Vuex from 'vuex';
import Vue from 'vue';
import todos from './modules/todos'
import users from './modules/users'
import blog from './modules/blog'

Vue.use(Vuex);

export default new Vuex.Store({
	modules: {
		todos,
		users,
		blog
	}
});
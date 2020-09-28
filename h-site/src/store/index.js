import Vue from "vue";
import Vuex from "vuex";
import users from "./modules/users";
import blog from "./modules/blog";
import projects from "./modules/projects"

Vue.use(Vuex);

export default new Vuex.Store({
  modules: {
    users,
    blog,
    projects,
  },
});

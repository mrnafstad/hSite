import db from "../../firebaseConfig.js";

const state = {
    projects: []
}

const getters = {
    projects: (state) => state.projects
}

const mutations = {
    newProject: (state, project) => state.projects.unshift(project)
}

const actions = {
    async getProjects({ commit }) {
        await db.projects.get().then(function(querySnapshot) {
            querySnapshot.forEach(function(doc) {
              //console.log(doc.id, "=>", doc.data())
              commit("newPost", {
                id: doc.id,
                name: doc.data().name,
                repoUrl: doc.data().repoUrl,
                siteUrl: doc.data().siteUrl,
                progress: doc.data().state,
              });
              console.log(doc.data().name)
            });
        })
    }
}

export default {
    state,
    getters,
    actions,
    mutations,
  };